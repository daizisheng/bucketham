% dualham.w -- literate; ctangle -> .c, cweave -> .pdf
\datethis
@* Introduction.
This program counts \&{open} knight's tours (open Hamiltonian paths) of the
$m\times n$ knight graph by a \&{dual\--frontier} transfer that is aggregated
with \&{bucket sorting} rather than a trie, so it parallelizes; and it extracts
the \&{periodic} transfer once the frontier stabilizes, so that far columns are
reached by a cheap sparse matrix\--vector product (SpMV) instead of re\--running
the generator.

The method is Knuth's (as in \.{DYNAHAM}/\.{dynahamp}), reimplemented from
scratch. Work in the augmented graph $G^+=G+K_1$: a new \&{apex} vertex joins
every cell, and an open tour of~$G$ is a Hamiltonian cycle of~$G^+$ through the
apex. Place the cells column by column. After $s$ cells are placed, the
\&{frontier} is the set of not\--yet\--placed vertices adjacent to a placed one,
plus the apex. Each frontier vertex is \&{bare} (degree~0), \&{outer} (degree~1,
a subpath endpoint), or \&{inner} (degree~2). The state is that pattern together
with the pairing of the outer endpoints. Because it is coded by \&{relative}
positions, not absolute vertex numbers, the code is translation\--invariant:
the sorted state set repeats with some small period once we are in the bulk,
which is exactly what makes the transfer periodic and the SpMV possible.

@c
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <sys/resource.h>
#include <unistd.h>
@#
@<Types@>@;
@<Globals@>@;
@<Subroutines@>@;
@<Main@>@;

@* Board and frontier.
Cell $(r,c)$ is vertex $v=cm+r$; there are $V=mn$ cells and the apex is~$V$.
|frontier_before(s)| lists (sorted) the frontier just before cell~$s$: the apex,
$s$ itself (an ``extended'' frontier, always present), and each future cell $>s$
with an already\--placed neighbour.

@<Types@>=
typedef unsigned long long u64;

@ @d MAXF 512
@d MAXLEV 40
@<Globals@>=
int m,n,V;
int NB[64*64][8], ND[64*64];
static const int KR[8]={-2,-2,-1,-1,1,1,2,2};
static const int KC[8]={-1,1,-2,2,-2,2,-1,1};
int fr[MAXF], ifrb[64*64+2];

@ All large allocations go through checked wrappers: on failure (e.g. hitting a
|ulimit -v| cap) they print and |exit| gracefully, so a run can never provoke
the kernel's OOM killer.

@<Subroutines@>=
void* xmalloc(size_t n){ void*p=malloc(n); if(!p&&n){ fprintf(stderr,"OOM: malloc %zu bytes failed\n",n); exit(3); } return p; }
void* xrealloc(void*q,size_t n){ void*p=realloc(q,n); if(!p&&n){ fprintf(stderr,"OOM: realloc %zu bytes failed\n",n); exit(3); } return p; }
void* xcalloc(size_t a,size_t b){ void*p=calloc(a,b); if(!p&&a&&b){ fprintf(stderr,"OOM: calloc %zu*%zu failed\n",a,b); exit(3); } return p; }

@ @<Subroutines@>=
void build_board(void){ int r,c,k; V=m*n;
  for(c=0;c<n;c++)for(r=0;r<m;r++){ int v=c*m+r; ND[v]=0;
    for(k=0;k<8;k++){ int rr=r+KR[k],cc=c+KC[k];
      if(rr>=0&&rr<m&&cc>=0&&cc<n) NB[v][ND[v]++]=cc*m+rr; } } }
int frontier_before(int s){ int q=0,v,u,k; static char inF[64*64+2];
  for(v=0;v<=V;v++) inF[v]=0; inF[V]=1; if(s<V) inF[s]=1;
  for(v=s+1;v<V;v++) for(k=0;k<ND[v];k++){ u=NB[v][k]; if(u<s){ inF[v]=1; break; } }
  for(v=0;v<=V;v++) if(inF[v]){ fr[q]=v; ifrb[v]=q; q++; }
  return q; }

@* Mate table, derived edges, key, completion.
The state is |mate[]| over frontier positions: |-2| bare, |-1| inner, else the
partner position (outer). |add_derived(i,j)| joins the subpaths ending at $i,j$;
if they are the two ends of one subpath a cycle forms (flagged in |cycle|, both
ends set inner). The key is one byte per position; equal (also column\--shifted)
patterns give equal keys. A cycle through the apex is a valid open $m'$\--tour
iff the covered (inner) cells are the contiguous prefix $\{0..m'-1\}$ and the
rest of the board frontier is bare; |completion_mp| returns that~$m'$.

@<Globals@>=
int mate[MAXF]; int cycle;
#pragma omp threadprivate(mate,cycle)

@ @<Subroutines@>=
int add_derived(int i,int j){ cycle=0;
  if(mate[i]==-1||mate[j]==-1) return 0;
  if(mate[i]==-2){ if(mate[j]==-2){ if(i==j){cycle=1;return 0;} mate[i]=j;mate[j]=i;return 1; }
    mate[i]=mate[j]; mate[mate[j]]=i; mate[j]=-1; return 1;
  } else if(mate[j]==-2){ mate[j]=mate[i]; mate[mate[i]]=j; mate[i]=-1; return 1; }
  else if(mate[i]!=j){ mate[mate[i]]=mate[j]; mate[mate[j]]=mate[i]; mate[i]=mate[j]=-1; return 1; }
  mate[i]=mate[j]=-1; cycle=1; return 0; }
int keyof(int q,unsigned char*key){ int i;
  for(i=0;i<q;i++) key[i]= mate[i]==-2?0 : mate[i]==-1?255 : (unsigned char)(1+mate[i]); return q; }
int completion_mp(int q,int apexpos,int*frn,int s){ int k,mp=s+1;
  if(mate[apexpos]!=-1) return 0;
  k=0; while(k<q && frn[k]!=V && mate[k]==-1 && frn[k]==mp){ mp++; k++; }
  for(;k<q;k++){ if(frn[k]==V) continue; if(mate[k]!=-2) return 0; } return mp; }

@* Bucket aggregation.
Each step emits successor records; we sort by key and merge equal keys. |src|
carries the emitting state's index, used only while recording the edge tables.

@<Types@>=
typedef struct { long off; int len; u64 w; long src; } Rec;

@ @<Globals@>=
unsigned char*kp; long kcap,kuse; Rec*rc; long nr,rcap; unsigned char*cb;
@#
#define NTMAX 128
unsigned char* tkp[NTMAX]; long tkuse[NTMAX], tkcap[NTMAX];
Rec* trc[NTMAX]; long tnr[NTMAX], trcap[NTMAX];
int rcmp(const void*A,const void*B){ const Rec*a=A,*b=B; int l=a->len<b->len?a->len:b->len;
  int d=memcmp(cb+a->off,cb+b->off,l); return d?d:a->len-b->len; }

@ @<Subroutines@>=
void emit(unsigned char*key,int len,u64 w,long src){
  int t=omp_get_thread_num();
  if(tkuse[t]+len>tkcap[t]){tkcap[t]=tkcap[t]*2+len+65536;tkp[t]=xrealloc(tkp[t],tkcap[t]);}
  if(tnr[t]>=trcap[t]){trcap[t]=trcap[t]*2+65536;trc[t]=xrealloc(trc[t],trcap[t]*sizeof(Rec));}
  memcpy(tkp[t]+tkuse[t],key,len); trc[t][tnr[t]].off=tkuse[t]; trc[t][tnr[t]].len=len;
  trc[t][tnr[t]].w=w; trc[t][tnr[t]].src=src; tkuse[t]+=len; tnr[t]++; }
void merge_pools(void){ int t; long r; nr=0; kuse=0;
  for(t=0;t<omp_get_max_threads();t++){ for(r=0;r<tnr[t];r++){ Rec*R=&trc[t][r];
    if(kuse+R->len>kcap){kcap=kcap*2+R->len+65536;kp=xrealloc(kp,kcap);}
    if(nr>=rcap){rcap=rcap*2+65536;rc=xrealloc(rc,rcap*sizeof(Rec));}
    memcpy(kp+kuse,tkp[t]+R->off,R->len);
    rc[nr].off=kuse; rc[nr].len=R->len; rc[nr].w=R->w; rc[nr].src=R->src; kuse+=R->len; nr++; }
    tnr[t]=0; tkuse[t]=0; } }
@#
int rcmp_r(const void*A,const void*B,void*arg){ unsigned char*base=arg; const Rec*a=A,*b=B;
  int l=a->len<b->len?a->len:b->len; int d=memcmp(base+a->off,base+b->off,l); return d?d:a->len-b->len; }
@#
/* sort each thread's records in parallel, then merge+reduce across them in one
   pass (keeping global sort order so periodic ids stay consistent), capturing
   the integer edge table when recording. */
long in_ncur_g;
/* Parallel splitter merge: sort each thread's records, sample P-1 splitter keys
   that cut the key space into P balanced ranges, binary-search each thread for
   the range boundaries, then merge+reduce each range INDEPENDENTLY in parallel,
   writing to a per-range output buffer. Equal keys never span a range (splitter
   lower_bound is consistent across threads), so the ranges just concatenate.
   Parallel AND memory-safe: no full sorted copy of the records is materialized. */
#define NRNG 256
int keycmp(unsigned char*a,int la,unsigned char*b,int lb){ int l=la<lb?la:lb; int d=memcmp(a,b,l); return d?d:la-lb; }
typedef struct{ unsigned char* k; int len; } Samp;
int sampcmp(const void*A,const void*B){ const Samp*a=A,*b=B; return keycmp(a->k,a->len,b->k,b->len); }
static unsigned char* rk[NRNG]; static long rkcap[NRNG], rkuse[NRNG];
static int* rkl[NRNG]; static u64* rw_[NRNG]; static long rcnt[NRNG], rgcap[NRNG];
static Edge* re[NRNG]; static long rne[NRNG], recap[NRNG];
static long bnd[NTMAX][NRNG+1];
void sort_reduce(int s,int rec){
  int NT=omp_get_max_threads(), t, r; long i;
#pragma omp parallel for schedule(dynamic,1)
  for(t=0;t<NT;t++) if(tnr[t]) qsort_r(trc[t],tnr[t],sizeof(Rec),rcmp_r,tkp[t]);
  long tot=0; for(t=0;t<NT;t++) tot+=tnr[t];
  if(getenv("SRDBG")){ int q5; fprintf(stderr,"trc[0] after qsort_r:"); for(q5=0;q5<10&&q5<tnr[0];q5++){ Rec*R=&trc[0][q5]; int kn=0,z; for(z=0;z<R->len;z++) kn=(kn<<8)|(tkp[0]+R->off)[z]; fprintf(stderr," %d",kn);} fprintf(stderr,"\n"); }
  int P=NT*8; if(P>NRNG) P=NRNG; if(P<1) P=1; if((long)P>tot && tot>0) P=(int)tot;
  int S=P*8; if((long)S>tot) S=(int)tot;
  static Samp* samp=0; static long sampcap=0; if(sampcap<S+1){ sampcap=S+1; samp=xrealloc(samp,sampcap*sizeof(Samp)); }
  long sc=0;
  for(t=0;t<NT && tot>0;t++){ long nt=tnr[t]; if(nt==0) continue; long take=(long)((double)S*nt/tot); if(take<1) take=1; long j;
    for(j=0;j<take && sc<S;j++){ long idx=(long)((double)j*nt/take); if(idx>=nt) idx=nt-1; Rec*R=&trc[t][idx]; samp[sc].k=tkp[t]+R->off; samp[sc].len=R->len; sc++; } }
  S=(int)sc; qsort(samp,S,sizeof(Samp),sampcmp);
  if(getenv("SRDBG")){ int q4; fprintf(stderr,"samp after qsort (S=%d):",S); for(q4=0;q4<S&&q4<16;q4++){ int kn=0,z; for(z=0;z<samp[q4].len;z++) kn=(kn<<8)|samp[q4].k[z]; fprintf(stderr," %d",kn);} fprintf(stderr,"\n"); }
  static Samp spl[NRNG]; int np=0, p;
  for(p=1;p<P;p++){ long si=(long)((double)p*S/P); if(si>=S) si=S-1; if(S>0) spl[np++]=samp[si]; }
  if(getenv("SRDBG")){ fprintf(stderr,"P=%d np=%d tot=%ld  splitters:",P,np,tot); int pp2; for(pp2=0;pp2<np&&pp2<12;pp2++){ int kn=0,z; for(z=0;z<spl[pp2].len;z++) kn=(kn<<8)|spl[pp2].k[z]; fprintf(stderr," %d",kn);} fprintf(stderr,"\n"); }
  for(t=0;t<NT;t++){ bnd[t][0]=0; bnd[t][P]=tnr[t];
    for(p=0;p<np;p++){ long lo=0,hi=tnr[t]; while(lo<hi){ long mid=(lo+hi)/2; Rec*R=&trc[t][mid];
      if(keycmp(tkp[t]+R->off,R->len,spl[p].k,spl[p].len)<0) lo=mid+1; else hi=mid; } bnd[t][p+1]=lo; } }
#pragma omp parallel for schedule(dynamic,1)
  for(r=0;r<P;r++){
    int hp[NTMAX]; long hpos[NTMAX]; int hn=0, tt;
    for(tt=0;tt<NT;tt++) hpos[tt]=bnd[tt][r];
#define RLESS(a,b) ({ Rec*Ra=&trc[a][hpos[a]],*Rb=&trc[b][hpos[b]]; keycmp(tkp[a]+Ra->off,Ra->len,tkp[b]+Rb->off,Rb->len)<0; })
    for(tt=0;tt<NT;tt++) if(hpos[tt]<bnd[tt][r+1]){ int c=hn++; hp[c]=tt;
      while(c>0){ int pp=(c-1)/2; if(RLESS(hp[c],hp[pp])){ int x=hp[c];hp[c]=hp[pp];hp[pp]=x; c=pp; } else break; } }
    rkuse[r]=0; rcnt[r]=0; rne[r]=0; int curlen=-1; long curpos=0; u64 sw=0;
    while(hn>0){ int mt=hp[0]; Rec*R=&trc[mt][hpos[mt]]; unsigned char*rkey=tkp[mt]+R->off; int rlen=R->len;
      int same=(curlen==rlen && memcmp(rk[r]+curpos,rkey,rlen)==0);
      if(!same){ if(curlen>=0){ if(rcnt[r]>=rgcap[r]){ rgcap[r]=rgcap[r]*2+1024; rkl[r]=xrealloc(rkl[r],rgcap[r]*sizeof(int)); rw_[r]=xrealloc(rw_[r],rgcap[r]*sizeof(u64)); } rkl[r][rcnt[r]]=curlen; rw_[r][rcnt[r]]=sw; rcnt[r]++; }
        if(rkuse[r]+rlen>rkcap[r]){ rkcap[r]=rkcap[r]*2+rlen+4096; rk[r]=xrealloc(rk[r],rkcap[r]); }
        curpos=rkuse[r]; memcpy(rk[r]+rkuse[r],rkey,rlen); rkuse[r]+=rlen; curlen=rlen; sw=0; }
      sw=red(sw+R->w);
      if(rec){ if(rne[r]>=recap[r]){ recap[r]=recap[r]*2+1024; re[r]=xrealloc(re[r],recap[r]*sizeof(Edge)); } re[r][rne[r]].src=(int)R->src; re[r][rne[r]].dst=(int)rcnt[r]; re[r][rne[r]].c=1; rne[r]++; }
      hpos[mt]++;
      if(hpos[mt]>=bnd[mt][r+1]) hp[0]=hp[--hn];
      { int c=0; while(1){ int l=2*c+1,r2=2*c+2,sm=c; if(l<hn&&RLESS(hp[l],hp[sm]))sm=l; if(r2<hn&&RLESS(hp[r2],hp[sm]))sm=r2; if(sm==c)break; int x=hp[c];hp[c]=hp[sm];hp[sm]=x; c=sm; } }
    }
    if(curlen>=0){ if(rcnt[r]>=rgcap[r]){ rgcap[r]=rgcap[r]*2+1024; rkl[r]=xrealloc(rkl[r],rgcap[r]*sizeof(int)); rw_[r]=xrealloc(rw_[r],rgcap[r]*sizeof(u64)); } rkl[r][rcnt[r]]=curlen; rw_[r][rcnt[r]]=sw; rcnt[r]++; }
  }
  if(getenv("SRDBG")){ int rr3; for(rr3=0;rr3<P&&rr3<12;rr3++){ int fk=-1,lk=-1,z; if(rcnt[rr3]>0){ fk=0; for(z=0;z<rkl[rr3][0];z++) fk=(fk<<8)|rk[rr3][z]; long lo2=0,q; for(q=0;q<rcnt[rr3]-1;q++) lo2+=rkl[rr3][q]; lk=0; for(z=0;z<rkl[rr3][rcnt[rr3]-1];z++) lk=(lk<<8)|rk[rr3][lo2+z]; } fprintf(stderr,"  range %d: rcnt=%ld first=%d last=%d\n",rr3,rcnt[rr3],fk,lk); } }
  static long base[NRNG+1], kbase[NRNG+1], ebase[NRNG+1];
  base[0]=kbase[0]=ebase[0]=0;
  for(r=0;r<P;r++){ base[r+1]=base[r]+rcnt[r]; kbase[r+1]=kbase[r]+rkuse[r]; ebase[r+1]=ebase[r]+rne[r]; }
  ncur=base[P];
  curoff=xrealloc(curoff,(ncur+1)*sizeof(long)); curkl=xrealloc(curkl,(ncur+1)*sizeof(int));
  curw=xrealloc(curw,(ncur+1)*sizeof(u64)); curkp=xrealloc(curkp,kbase[P]+1);
#pragma omp parallel for schedule(dynamic,1)
  for(r=0;r<P;r++){ long off=kbase[r], j; memcpy(curkp+kbase[r], rk[r], rkuse[r]);
    for(j=0;j<rcnt[r];j++){ curkl[base[r]+j]=rkl[r][j]; curw[base[r]+j]=rw_[r][j]; curoff[base[r]+j]=off; off+=rkl[r][j]; } }
  for(t=0;t<NT;t++){ tnr[t]=0; tkuse[t]=0; }
  if(rec){ static long rne2[NRNG], eb2[NRNG+1]; int r2;
    /* coalesce each range's edges IN PLACE (parallel), giving global dsts; ranges
       hold disjoint, ordered dst ranges, so concatenating them is already globally
       (dst,src)-sorted -- no big duplicate array, no serial global sort. */
#pragma omp parallel for schedule(dynamic,1)
    for(r2=0;r2<P;r2++){ long j; for(j=0;j<rne[r2];j++) re[r2][j].dst += base[r2];
      qsort(re[r2],rne[r2],sizeof(Edge),ecmp);
      long o=0,pp; for(pp=0;pp<rne[r2];){ long q2=pp+1; u64 cc=1; while(q2<rne[r2]&&re[r2][q2].src==re[r2][pp].src&&re[r2][q2].dst==re[r2][pp].dst){cc++;q2++;} re[r2][o]=re[r2][pp]; re[r2][o].c=cc; o++; pp=q2; } rne2[r2]=o; }
    long tote=0; for(r2=0;r2<P;r2++) tote+=rne2[r2];
    eb2[0]=0; for(r2=0;r2<P;r2++) eb2[r2+1]=eb2[r2]+rne2[r2];
    Edge* eb=xmalloc((tote+1)*sizeof(Edge));
#pragma omp parallel for schedule(dynamic,1)
    for(r2=0;r2<P;r2++) memcpy(eb+eb2[r2], re[r2], rne2[r2]*sizeof(Edge));
    edges[reclev]=eb; nedge[reclev]=tote;
    nstate[reclev]=in_ncur_g; nstate[reclev+1]=ncur;
    comps[reclev]=xmalloc((reccomp_n+1)*sizeof(Comp)); memcpy(comps[reclev],reccomp_buf,reccomp_n*sizeof(Comp)); ncomp[reclev]=reccomp_n;
    reclev++; }
}

@* The transfer step.
The bucket holds states as (key,weight) over the previous frontier. To place
cell~$s$: build the base mate table |bmate| over the new frontier (survivors
carried, newcomers bare) with $s$ in a temporary slot, then bring $s$ to
degree~2 by adding its $2-\deg(s)$ edges to future neighbours or the apex. When
|rec| is set we also record the integer edge table (src index $\to$ dst index)
and the completion contributions for this level.

@<Globals@>=
unsigned char*curkp; long*curoff; int*curkl; u64*curw; long ncur;
u64 cnt[1<<16];
u64 MODP=0;   /* 0 = exact u64; set to a prime for one CRT residue */
static inline u64 red(u64 x){ return MODP? x%MODP : x; }
int qnew, posS, apexnew, STEMP; int bmate[MAXF], o2n[MAXF];
#pragma omp threadprivate(bmate,STEMP)
@#
int recording, reclev;
typedef struct{ int src,dst; u64 c; } Edge;
typedef struct{ int src,delta; u64 mult; } Comp;
Edge* edges[MAXLEV]; long nedge[MAXLEV];
Comp* comps[MAXLEV]; long ncomp[MAXLEV];
long nstate[MAXLEV+1]; int Plevs;
int period_g, c0_g, recend_col_g;   /* saved for dumping the extracted tables */
int stop_after_record=0; int direct_cov_g=0;
static u64 colfp_g[4096]; const char* ckpt_path=0; int ckpt_every=4;   /* checkpoint every K columns */
u64* seedv; long seedn;
Comp* reccomp_buf; long reccomp_n, reccomp_cap;
int ecmp(const void*A,const void*B){ const Edge*a=A,*b=B; if(a->dst!=b->dst) return a->dst-b->dst; return a->src-b->src; }

@ @<Subroutines@>=
void build_bmate(int*omate,int qold){ int i;
  for(i=0;i<=qnew;i++) bmate[i]=-2; STEMP=qnew;
  for(i=0;i<qold;i++){ int dst=(i==posS)?STEMP:o2n[i]; if(dst<0) continue;
    if(omate[i]==-1) bmate[dst]=-1;
    else if(omate[i]>=0){ int op=omate[i]; int pdst=(op==posS)?STEMP:o2n[op]; bmate[dst]=pdst>=0?pdst:-2; } } }

@ @<Subroutines@>=
void run_step(int s,int rec){
  int i; long si; long in_ncur=ncur; int qold=frontier_before(s); posS=ifrb[s];
  static int frold[MAXF]; for(i=0;i<qold;i++) frold[i]=fr[i];
  qnew=frontier_before(s+1); static int ifrnew[64*64+2], frnew[MAXF];
  for(i=0;i<qnew;i++) frnew[i]=fr[i]; for(i=0;i<=V;i++) ifrnew[i]=-1; for(i=0;i<qnew;i++) ifrnew[frnew[i]]=i;
  apexnew=ifrnew[V];
  int nbr[16],rr=0; nbr[rr++]=apexnew;
  for(i=0;i<ND[s];i++){ int w=NB[s][i]; if(w>s && ifrnew[w]>=0) nbr[rr++]=ifrnew[w]; }
  for(i=0;i<qold;i++) o2n[i]=(frold[i]==s)?-1:ifrnew[frold[i]];
  in_ncur_g=ncur; reccomp_n=0;
  @<Expand every state at cell |s|@>;
  @<Sort, reduce, and (if recording) capture the edge table@>;
}

@ Each state expands independently, so the loop runs in parallel (except while
recording, when it stays serial to keep the completion log ordered). The
transition scratch (|mate|, |bmate|, |cycle|, |STEMP|) is |threadprivate|; each
thread emits into its own pool; |cnt| is updated in a critical section.

@<Expand every state at cell |s|@>=
#pragma omp parallel for schedule(dynamic,32)
for(si=0;si<ncur;si++){
  int i,a,b,nl,deg,need,mp; unsigned char*ok=curkp+curoff[si]; int okl=curkl[si]; u64 w=curw[si];
  int omate[MAXF]; unsigned char nk[MAXF];
  for(i=0;i<okl;i++){int cc=ok[i]; omate[i]= cc==0?-2 : cc==255?-1 : cc-1;}
  deg=omate[posS]==-2?0:omate[posS]==-1?2:1; need=2-deg;
  if(need==0){ build_bmate(omate,qold); for(i=0;i<qnew;i++) mate[i]=bmate[i]; nl=keyof(qnew,nk); emit(nk,nl,w,si); }
  else if(need==1){ for(a=0;a<rr;a++){ build_bmate(omate,qold); for(i=0;i<=qnew;i++) mate[i]=bmate[i];
    if(add_derived(STEMP,nbr[a])){ nl=keyof(qnew,nk); emit(nk,nl,w,si); }
    else if(cycle){ mp=completion_mp(qnew,apexnew,frnew,s); if(mp) @<Credit completion@>; } } }
  else { for(a=0;a<rr;a++)for(b=a+1;b<rr;b++){ build_bmate(omate,qold); for(i=0;i<=qnew;i++) mate[i]=bmate[i];
    if(!add_derived(STEMP,nbr[a])){ if(cycle){ mp=completion_mp(qnew,apexnew,frnew,s); if(mp) @<Credit completion@>; } continue; }
    if(add_derived(STEMP,nbr[b])){ nl=keyof(qnew,nk); emit(nk,nl,w,si); }
    else if(cycle){ mp=completion_mp(qnew,apexnew,frnew,s); if(mp) @<Credit completion@>; } } }
}

@ @<Credit completion@>=
{
#pragma omp critical
  { cnt[mp]=red(cnt[mp]+w);
    if(rec){ if(reccomp_n>=reccomp_cap){ reccomp_cap=reccomp_cap*2+(1<<20); reccomp_buf=xrealloc(reccomp_buf,reccomp_cap*sizeof(Comp)); }
      reccomp_buf[reccomp_n].src=si; reccomp_buf[reccomp_n].delta=mp-(s+1); reccomp_buf[reccomp_n].mult=1; reccomp_n++; } }
}

@ The sort and reduce run across the per-thread pools in parallel (see
|sort_reduce|); |in_ncur_g| records the input level size for the edge table.

@<Sort, reduce, and (if recording) capture the edge table@>=
sort_reduce(s,rec);

@* Checkpointing the build.
The build (the sweep that discovers the periodic transfer) is the long part; we
checkpoint it at every column boundary before recording begins, dumping the
current bucket, the running counts, and the fingerprints. A build that dies
resumes from the last boundary; the recording columns themselves are few, so if
a crash lands inside them we simply restart from the last pre\--recording
checkpoint.

@<Subroutines@>=
void save_ckpt(int nexts){ if(!ckpt_path) return; FILE*f=fopen(ckpt_path,"wb");
  fwrite(&m,sizeof(int),1,f); fwrite(&nexts,sizeof(int),1,f);
  fwrite(cnt,sizeof(u64),1<<16,f); fwrite(colfp_g,sizeof(u64),4096,f);
  fwrite(&ncur,sizeof(long),1,f); fwrite(&kuse,sizeof(long),1,f);
  fwrite(curoff,sizeof(long),ncur,f); fwrite(curkl,sizeof(int),ncur,f);
  fwrite(curw,sizeof(u64),ncur,f);
  { long kb=0,i; for(i=0;i<ncur;i++) kb+=curkl[i]; fwrite(&kb,sizeof(long),1,f); fwrite(curkp,1,kb,f); }
  fclose(f); }
int load_ckpt(const char*path){ FILE*f=fopen(path,"rb"); if(!f) return -1; int nexts,mm;
  if(fread(&mm,sizeof(int),1,f)!=1) return -1; m=mm; fread(&nexts,sizeof(int),1,f);
  fread(cnt,sizeof(u64),1<<16,f); fread(colfp_g,sizeof(u64),4096,f);
  fread(&ncur,sizeof(long),1,f); fread(&kuse,sizeof(long),1,f);
  curoff=xrealloc(curoff,(ncur+1)*sizeof(long)); curkl=xrealloc(curkl,(ncur+1)*sizeof(int)); curw=xrealloc(curw,(ncur+1)*sizeof(u64));
  fread(curoff,sizeof(long),ncur,f); fread(curkl,sizeof(int),ncur,f); fread(curw,sizeof(u64),ncur,f);
  { long kb; fread(&kb,sizeof(long),1,f); curkp=xrealloc(curkp,kb+1); fread(curkp,1,kb,f); }
  fclose(f); fprintf(stderr,"resumed build from %s at cell %d (col %d)\n",path,nexts,nexts/mm); return nexts; }

@ @* Dumping and reloading the extracted tables.
The expensive part is the build that discovers and records the periodic transfer.
Once recorded, the tables (edge tables, completions, seed vector, level sizes,
period, and the directly\--counted columns up to |recend_col_g|) are small and
self\--contained, so we dump them to disk. A later run reloads them and does only
the cheap SpMV --- and a build that dies can be re\--run without touching the SpMV.

@<Subroutines@>=
void dump_tables(const char*path){ FILE*f=fopen(path,"wb"); int L;
  int hdr[5]={m,period_g,c0_g,Plevs,recend_col_g}; fwrite(hdr,sizeof(int),5,f);
  fwrite(&seedn,sizeof(long),1,f); fwrite(seedv,sizeof(u64),seedn,f);
  fwrite(nstate,sizeof(long),Plevs+1,f);
  for(L=0;L<Plevs;L++){ fwrite(&nedge[L],sizeof(long),1,f); fwrite(edges[L],sizeof(Edge),nedge[L],f);
    fwrite(&ncomp[L],sizeof(long),1,f); fwrite(comps[L],sizeof(Comp),ncomp[L],f); }
  int nc=(recend_col_g+1)*m; fwrite(&nc,sizeof(int),1,f); fwrite(cnt,sizeof(u64),nc,f);
  fclose(f); fprintf(stderr,"dumped tables to %s (period %d, c0 %d, %d levels)\n",path,period_g,c0_g,Plevs); }
int load_tables(const char*path){ FILE*f=fopen(path,"rb"); if(!f) return 0; int L,hdr[5];
  if(fread(hdr,sizeof(int),5,f)!=5) return 0; m=hdr[0]; period_g=hdr[1]; c0_g=hdr[2]; Plevs=hdr[3]; recend_col_g=hdr[4];
  fread(&seedn,sizeof(long),1,f); seedv=xmalloc(seedn*sizeof(u64)); fread(seedv,sizeof(u64),seedn,f);
  fread(nstate,sizeof(long),Plevs+1,f);
  for(L=0;L<Plevs;L++){ fread(&nedge[L],sizeof(long),1,f); edges[L]=xmalloc((nedge[L]+1)*sizeof(Edge)); fread(edges[L],sizeof(Edge),nedge[L],f);
    fread(&ncomp[L],sizeof(long),1,f); comps[L]=xmalloc((ncomp[L]+1)*sizeof(Comp)); fread(comps[L],sizeof(Comp),ncomp[L],f); }
  int nc; fread(&nc,sizeof(int),1,f); memset(cnt,0,sizeof(cnt)); fread(cnt,sizeof(u64),nc,f);
  fclose(f); return 1; }

@ @* Driver: build, extract the period, then SpMV.
Sweep an $m\times W_b$ strip; when the boundary key\--set (fingerprinted) repeats
with period $p\in\{1,2\}$ at column~$c_0$, record the next $p$ columns' edge
tables and capture the seed vector. Then iterate the SpMV to reach column
$N$: at each recorded level apply the edge table to the state vector and add the
level's completions to |cnt2|. The two agree on every column the SpMV fully
covers ($c\ge c_0+3$), and the SpMV alone yields the columns past~$W_b$.

@<Subroutines@>=
int resume_from=0;
void seed_bucket(void){ int i; int q=frontier_before(0);
  for(i=0;i<q;i++) mate[i]=-2; unsigned char key[MAXF]; int kl=keyof(q,key);
  curkp=xmalloc(1<<20); curoff=xmalloc(sizeof(long)); curkl=xmalloc(sizeof(int)); curw=xmalloc(sizeof(u64));
  memcpy(curkp,key,kl); curoff[0]=0; curkl[0]=kl; curw[0]=1; ncur=1; }
@#
u64 cnt2[1<<16];
int run_periodic(int mm,int Wb,int Next){
  int i,s,c; m=mm; n=Wb; build_board(); memset(cnt2,0,sizeof(cnt2));
  recording=0; reclev=0;
  if(resume_from<=0){ memset(cnt,0,sizeof(cnt)); seed_bucket(); }  /* fresh; else state is loaded */
  int c0=-1,period=0,recstart=-1,recend=-1;
  for(s=(resume_from>0?resume_from:0);s<V;s++){
    if(recording && s==recstart){ seedn=ncur; seedv=xmalloc(ncur*sizeof(u64)); for(i=0;i<ncur;i++) seedv[i]=curw[i]; nstate[0]=ncur; reclev=0; }
    run_step(s, recording && s>=recstart && s<=recend);
    if(s%m==m-1){ c=s/m; u64 h=1469598103934665603ULL; long t;
      for(t=0;t<ncur;t++){ unsigned char*kk=curkp+curoff[t]; int L=curkl[t],z; for(z=0;z<L;z++){ h^=kk[z]; h*=1099511628211ULL; } h^=0x9e; h*=1099511628211ULL; }
      colfp_g[c]=h;
      if(c0<0){ if(c>=1&&colfp_g[c]==colfp_g[c-1]){period=1;c0=c;} else if(c>=2&&colfp_g[c]==colfp_g[c-2]){period=2;c0=c;}
        if(c0>=0){ recstart=(c0+1)*m; recend=(c0+1+period)*m-1; Plevs=period*m; recording=1;
          period_g=period; c0_g=c0; recend_col_g=c0+period;
          fprintf(stderr,"stable col %d, period %d\n",c0,period); } } }
    if(s%m==m-1 && getenv("DBG")) fprintf(stderr,"  col %d ncur=%ld rec=%d\n",s/m,ncur,recording);
    if(!recording && s%m==m-1 && (s/m)%ckpt_every==0) save_ckpt(s+1);
    if(recording && s==recend && stop_after_record) break;
  }
  direct_cov_g = stop_after_record ? recend_col_g : n;
  spmv_run(Next, 1);
  return c0;
}
@#
/* SpMV from the (built or loaded) tables out to column |Nto|. |crosscheck|
   compares against the direct |cnt| where both are valid. Uses |cnt2| scratch. */
int build_extract(int mm,int Wb){ return run_periodic(mm,Wb,Wb); }
void spmv_run(int Nto,int crosscheck){
  int i; memset(cnt2,0,sizeof(cnt2));
  if(MODP){ long j; for(j=0;j<seedn;j++) seedv[j]=red(seedv[j]); }
  @<Iterate the SpMV out to column $N$@>;
  { int good=1,cc;
    if(crosscheck) for(cc=c0_g+3;cc<=Nto && cc<=direct_cov_g;cc++) if(cnt[cc*m]!=cnt2[cc*m]){ good=0;
      fprintf(stderr,"MISMATCH c=%d direct=%llu spmv=%llu\n",cc,cnt[cc*m],cnt2[cc*m]); }
    for(cc=1;cc<=Nto;cc++){ u64 val = cc<=recend_col_g? red(cnt[cc*m]) : cnt2[cc*m]; if(val) printf("open %dx%d = %llu%s\n",m,cc,val, cc>recend_col_g?"  (SpMV)":""); }
    if(crosscheck) fprintf(stderr,"%s\n", good?"SpMV matches direct":"SpMV MISMATCH"); }
}

@ @<Iterate the SpMV out to column $N$@>=
{ u64* v=xmalloc(nstate[0]*sizeof(u64)); memcpy(v,seedv,nstate[0]*sizeof(u64));
  int basecol=c0_g+1;
  while(basecol<=Nto){ int L;
    for(L=0;L<Plevs;L++){ int abscol=basecol+L/m, substep=L%m; long e;
      for(e=0;e<ncomp[L];e++){ int idx=abscol*m+substep+1+comps[L][e].delta; cnt2[idx]=red(cnt2[idx]+v[comps[L][e].src]*comps[L][e].mult); }
      u64* vn=xcalloc(nstate[L+1],sizeof(u64));
      { long ne=nedge[L]; int NT=omp_get_max_threads(),tt; static long spl[NTMAX+1];
        spl[0]=0; spl[NT]=ne;
        for(tt=1;tt<NT;tt++){ long g=(long)((double)tt*ne/NT);
          while(g>0&&g<ne&&edges[L][g].dst==edges[L][g-1].dst) g++; spl[tt]=g; }
#pragma omp parallel for schedule(dynamic,1)
        for(tt=0;tt<NT;tt++){ long e2; for(e2=spl[tt];e2<spl[tt+1];e2++){ int d=edges[L][e2].dst;
          vn[d]=red(vn[d]+v[edges[L][e2].src]*edges[L][e2].c); } } }
      free(v); v=vn; }
    basecol+=period_g; }
  free(v); }

@ @<Report and cross\--check@>=
{ int good=1,cc;
  for(cc=c0+3;cc<=Next && cc<=Wb;cc++) if(cnt[cc*m]!=cnt2[cc*m]){ good=0;
    fprintf(stderr,"MISMATCH c=%d direct=%llu spmv=%llu\n",cc,cnt[cc*m],cnt2[cc*m]); }
  for(cc=1;cc<=Next;cc++){ u64 val = cc<=Wb? cnt[cc*m] : cnt2[cc*m]; if(val) printf("open %dx%d = %llu%s\n",m,cc,val, cc>Wb?"  (SpMV)":""); }
  fprintf(stderr,"%s\n", good?"SpMV matches direct":"SpMV MISMATCH"); }

@ A fast standalone correctness harness for |sort_reduce|: fill the thread pools
with |N| synthetic records whose keys encode one of |K| numbers in |L| big-endian
bytes (so lexicographic == numeric) with random small weights, run the merge, and
compare the reduced (key,weight) output against a direct per-number weight sum.
Exercises many duplicate keys spread across threads and splitter ranges.

@<Subroutines@>=
void tpush(int t,unsigned char*key,int len,u64 w){
  if(tkuse[t]+len>tkcap[t]){tkcap[t]=tkcap[t]*2+len+65536;tkp[t]=xrealloc(tkp[t],tkcap[t]);}
  if(tnr[t]>=trcap[t]){trcap[t]=trcap[t]*2+65536;trc[t]=xrealloc(trc[t],trcap[t]*sizeof(Rec));}
  memcpy(tkp[t]+tkuse[t],key,len); trc[t][tnr[t]].off=tkuse[t]; trc[t][tnr[t]].len=len;
  trc[t][tnr[t]].w=w; trc[t][tnr[t]].src=(long)tnr[t]; tkuse[t]+=len; tnr[t]++; }
int test_merge(long N,int L,int K,int rec){
  int NT=omp_get_max_threads(), t; long j; unsigned long seed=88172645463325252UL;
#define RND (seed^=seed<<13, seed^=seed>>7, seed^=seed<<17, seed)
  u64* ref=xcalloc(K,sizeof(u64)); MODP=0;
  for(j=0;j<N;j++){ unsigned long r=RND; int kn=(int)(r%K); u64 w=1+(RND%9);
    unsigned char key[64]; int i; for(i=0;i<L;i++) key[i]=((unsigned long)kn>>(8*(L-1-i)))&0xff;
    int t2=(int)(RND%NT); tpush(t2,key,L,w); ref[kn]=ref[kn]+w; }
  reccomp_n=0; recording=1; reclev=0; nstate[0]=0; sort_reduce(0,rec);
  /* compare */
  int bad=0; long i; long prev=-1;
  for(i=0;i<ncur;i++){ unsigned char*k=curkp+curoff[i]; long kn=0,z; for(z=0;z<curkl[i];z++) kn=(kn<<8)|k[z];
    if(curkl[i]!=L){ printf("len mismatch\n"); bad++; break; }
    if(kn<=prev){ printf("NOT SORTED at %ld: kn=%ld prev=%ld\n",i,kn,prev); bad++; break; } prev=kn;
    if((u64)curw[i]!=ref[kn]){ printf("WEIGHT mismatch kn=%ld got=%llu ref=%llu\n",kn,(u64)curw[i],ref[kn]); bad++; if(bad>5)break; } }
  long present=0; for(j=0;j<K;j++) if(ref[j]) present++;
  if(ncur!=present){ printf("COUNT mismatch: ncur=%ld distinct-present=%ld\n",ncur,present); bad++; }
  printf("test N=%ld L=%d K=%d rec=%d NT=%d: %s (ncur=%ld)\n",N,L,K,rec,NT,bad?"FAIL":"OK",ncur);
  free(ref); return bad; }

@ 
@* Main.
No arguments: self\--checks (direct sweep vs known values, and SpMV vs direct).
|dualham m Wb N|: build to column $W_b$, then extend to column~$N$ by SpMV.

@<Main@>=
void spmv_run(int Nto,int crosscheck);
int build_extract(int mm,int Wb);
int test_merge(long,int,int,int);
int main(int argc,char*argv[]){
  if(argc>=2 && !strcmp(argv[1],"test")){ return test_merge(argc>2?atol(argv[2]):200000, argc>3?atoi(argv[3]):3, argc>4?atoi(argv[4]):700, argc>5?atoi(argv[5]):0); }
  { long ram=sysconf(_SC_PHYS_PAGES)*(long)sysconf(_SC_PAGE_SIZE); long cap=(long)(ram*0.85);
    char*e=getenv("MEMCAP_GB"); if(e) cap=atol(e)*(1L<<30);
    struct rlimit rl={cap,cap}; setrlimit(RLIMIT_AS,&rl);
    fprintf(stderr,"self memory cap: %.0f GB (set MEMCAP_GB to override)\n",cap/1073741824.0); }
  if(argc>=5 && !strcmp(argv[1],"build")){ /* build m Wb file [ckpt] : build+extract, dump tables */
    stop_after_record=1; if(argc>=6) ckpt_path=argv[5]; if(argc>=7) ckpt_every=atoi(argv[6]);
    build_extract(atoi(argv[2]),atoi(argv[3])); dump_tables(argv[4]); return 0; }
  if(argc==5 && !strcmp(argv[1],"resume")){ /* resume ckpt Wb file : continue a build */
    stop_after_record=1; resume_from=load_ckpt(argv[2]); ckpt_path=argv[2];
    build_extract(m,atoi(argv[3])); dump_tables(argv[4]); return 0; }
  if(argc>=4 && !strcmp(argv[1],"run")){ /* run file Nto [prime] : load tables, SpMV (mod prime if given) */
    if(!load_tables(argv[2])){ fprintf(stderr,"cannot load %s\n",argv[2]); return 1; }
    if(argc>=5) MODP=strtoull(argv[4],0,10);
    spmv_run(atoi(argv[3]),0); return 0; }
  if(argc==4){ run_periodic(atoi(argv[1]),atoi(argv[2]),atoi(argv[3])); return 0; }
  struct{int m,Wb,c;u64 e;} chk[]={{5,12,4,82},{5,12,6,18784},{5,12,8,18061054ULL},
    {5,12,9,264895640ULL},{5,12,10,7886117822ULL},{7,6,4,6378},{6,7,6,3318960},{0,0,0,0}};
  int i,bad=0,lm=0,lw=0;
  for(i=0;chk[i].m;i++){
    if(chk[i].m!=lm||chk[i].Wb!=lw){ run_periodic(chk[i].m,chk[i].Wb,chk[i].Wb); lm=chk[i].m; lw=chk[i].Wb; }
    u64 g=cnt[chk[i].c*chk[i].m], e=red(chk[i].e);
    printf("open %dx%d = %llu  exp(mod p) %llu  %s\n",chk[i].m,chk[i].c,g,e,g==e?"OK":"FAIL");
    if(g!=e) bad++; }
  printf("%s\n",bad?"SOME FAILED":"ALL OK"); return bad?1:0;
}
