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
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <sys/resource.h>
#include <unistd.h>
#include <pthread.h>
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
typedef unsigned __int128 u128;   /* exact build weights (u64 seeds overflow at m>=7) */
typedef unsigned int u32;

@ @d MAXF 512
@d MAXLEV 400
@<Globals@>=
int m,n,V;
int NB[64*64][8], ND[64*64];
static const int KR[8]={-2,-2,-1,-1,1,1,2,2};
static const int KC[8]={-1,1,-2,2,-2,2,-1,1};
int fr[MAXF], ifrb[64*64+2];

@ All large allocations go through checked wrappers: on failure they print and
|exit| gracefully, so a run can never provoke the kernel's OOM killer.

@<Subroutines@>=
void* xmalloc(size_t n){ void*p=malloc(n); if(!p&&n){ fprintf(stderr,"OOM: malloc %zu bytes failed\n",n); exit(3); } return p; }
void* xrealloc(void*q,size_t n){ void*p=realloc(q,n); if(!p&&n){ fprintf(stderr,"OOM: realloc %zu bytes failed\n",n); exit(3); } return p; }
void* xcalloc(size_t a,size_t b){ void*p=calloc(a,b); if(!p&&a&&b){ fprintf(stderr,"OOM: calloc %zu*%zu failed\n",a,b); exit(3); } return p; }

@ To bound {\it physical\/} memory -- the resident set size, which is what
actually triggers the kernel OOM killer -- a detached watcher thread samples the
RSS and |_exit(3)|s gracefully if it exceeds the cap. We deliberately do {\bf not}
cap the virtual address space with |RLIMIT_AS|: glibc arenas, per\--thread stacks,
and |qsort|'s scratch |mmap| reserve far more virtual than resident memory, so an
|RLIMIT_AS| cap makes those (non-|xmalloc|) allocations fail with a hard crash
instead of a clean exit. RSS is the quantity that matters for OOM, and the largest
single allocation here (an edge table, a few GB) fits comfortably inside the
headroom between the 85\%-RAM cap and total RAM, so the 0.5\,s sampling never lets
resident memory overshoot the machine.

@<Globals@>=
long mem_cap_bytes;

@ @<Subroutines@>=
static long rss_bytes(void){ FILE*f=fopen("/proc/self/statm","r"); if(!f) return 0;
  long sz=0,res=0; if(fscanf(f,"%ld %ld",&sz,&res)!=2) res=0; fclose(f);
  return res*(long)sysconf(_SC_PAGE_SIZE); }
static void* mem_watcher(void*arg){ (void)arg;
  for(;;){ long r=rss_bytes();
    if(mem_cap_bytes>0 && r>mem_cap_bytes){
      fprintf(stderr,"\nMEMORY CAP: RSS %.1f GB exceeds cap %.1f GB -- exiting"
        " gracefully (raise MEMCAP_GB to allow more)\n",
        r/1073741824.0, mem_cap_bytes/1073741824.0); fflush(stderr); _exit(3); }
    usleep(300000); }
  return 0; }
/* The binding limit is the tightest of physical RAM and the cgroup |memory.max|
   over our cgroup-v2 hierarchy: a container/slice cap is often well below
   physical RAM, and the kernel's cgroup OOM killer fires at THAT limit, not the
   physical one -- so a watcher that only knew physical RAM would be OOM-killed
   before it ever tripped. */
static long cgroup_mem_limit(void){
  char line[4096],*p; char path[4096]={0};
  FILE*f=fopen("/proc/self/cgroup","r"); if(!f) return -1;
  while(fgets(line,sizeof line,f)) if((p=strstr(line,"0::"))){ p+=3;
    char*nl=strchr(p,'\n'); if(nl)*nl=0; snprintf(path,sizeof path,"%s",p); break; }
  fclose(f); if(!path[0]) return -1;
  long lim=-1; char base[4200]; snprintf(base,sizeof base,"/sys/fs/cgroup%s",path);
  int up; for(up=0; up<16; up++){
    char mm[4300]; snprintf(mm,sizeof mm,"%s/memory.max",base);
    FILE*g=fopen(mm,"r");
    if(g){ char v[64]={0}; if(fgets(v,sizeof v,g) && strncmp(v,"max",3)){
      long x=atol(v); if(x>0 && (lim<0||x<lim)) lim=x; } fclose(g); }
    if(!strcmp(base,"/sys/fs/cgroup")) break;
    char*s=strrchr(base,'/'); if(!s) break; *s=0;
    if(strlen(base) < strlen("/sys/fs/cgroup")) break; }
  return lim; }
void start_mem_watcher(void){
  long ram=sysconf(_SC_PHYS_PAGES)*(long)sysconf(_SC_PAGE_SIZE);
  long cg=cgroup_mem_limit();
  long avail = (cg>0 && cg<ram) ? cg : ram;   /* tightest binding limit */
  mem_cap_bytes=(long)(avail*0.85);
  char*e=getenv("MEMCAP_GB"); if(e) mem_cap_bytes=atol(e)*(1L<<30);
  pthread_t th; pthread_create(&th,0,mem_watcher,0); pthread_detach(th);
  fprintf(stderr,"self memory cap: %.0f GB RSS (%s%s; set MEMCAP_GB to override)\n",
    mem_cap_bytes/1073741824.0,
    (cg>0&&cg<ram)?"cgroup-limited":"85% of physical",
    (cg>0&&cg<ram)?"" : ""); }

@ The sweep order is the vertex numbering: cell processed at column~$c$,
in\--column position~$p$ is vertex $v=cm+p$. |rowat[p]| is the actual board row
at position~$p$ and |posof[]| its inverse. By default position~$=$~row
(top\--to\--bottom). With |OUTIN| set we number positions {\it outside\--in\/}
($0,m-1,1,m-2,\ldots$) so that after each completed symmetric pair the
processed\--within\--column row set is closed under the reflection $r\to m-1-r$
--- the precondition for the reflection quotient to be applied at more than just
the column boundary.

@<Globals@>=
int rowat[64], posof[64];

@ @<Subroutines@>=
void build_board(void){ int c,k,p; V=m*n;
  if(getenv("OUTIN")){ int lo=0,hi=m-1,i; for(p=0;p<m;p++){ if(p%2==0) rowat[p]=lo++; else rowat[p]=hi--; } (void)i; }
  else { for(p=0;p<m;p++) rowat[p]=p; }
  for(p=0;p<m;p++) posof[rowat[p]]=p;
  { int pp,rr; for(p=0;p<m;p++){ int ok=1; static char inset[64];
      for(rr=0;rr<m;rr++) inset[rr]=0; for(pp=0;pp<=p;pp++) inset[rowat[pp]]=1;
      for(rr=0;rr<m;rr++) if(inset[rr] && !inset[m-1-rr]){ ok=0; break; } sym_at[p]=ok; } }
  if(getenv("DBGROW")){ fprintf(stderr,"rowat="); for(p=0;p<m;p++) fprintf(stderr,"%d ",rowat[p]); fprintf(stderr,"  sym_at="); for(p=0;p<m;p++) fprintf(stderr,"%d",sym_at[p]); fprintf(stderr,"\n"); }
  for(c=0;c<n;c++)for(p=0;p<m;p++){ int r0=rowat[p]; int v=c*m+p; ND[v]=0;
    for(k=0;k<8;k++){ int rr=r0+KR[k],cc=c+KC[k];
      if(rr>=0&&rr<m&&cc>=0&&cc<n) NB[v][ND[v]++]=cc*m+posof[rr]; } } }
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
/* The record carries its key INLINE (keys are $\le 2m+1$ bytes; |KEYMAX| holds
   $m\le 9$), so there is no separate key pool to allocate/copy and the sort's
   |memcmp| hits contiguous memory instead of chasing a pointer. The weight is
   split into two |u64| halves so the struct is 8-byte, not 16-byte, aligned --
   a u128 member would pad the record to 48 bytes. All records of one level share
   the same key length |g_keylen| (=|qnew|), so no per-record |len| is stored.
   Net: 48 bytes + a 13-byte pooled key -> 40 bytes self-contained. */
#define KEYMAX 20
/* Weight-free: the record carries only its emitting state's index (for the edge
   table) and its key. No weight -- the build extracts the edge/completion
   structure only, and the SpMV replay from col 0 regenerates every weight.
   24 bytes, and the state set drops its weight column too. */
typedef struct { int src; unsigned char key[KEYMAX]; } Rec;

@ @<Globals@>=
long kuse;   /* vestigial: still written to the checkpoint stream (always 0) */
@#
#define NTMAX 128
Rec* trc[NTMAX]; long tnr[NTMAX], trcap[NTMAX];
int g_keylen;   /* key length (=qnew) for the level currently being reduced */

@ Hash aggregation (experimental, env |HASHAGG|; serial for now). Instead of
emitting a full record per successor and sorting, each successor key is hashed on
the fly to a distinct\--state id and only an 8\--byte |(dst,src)| edge is kept --
so the emit pool (the build's memory peak) becomes an 8B/edge accumulator instead
of a 24B/record one. The distinct keys are sorted once at the end so state ids
stay canonical (the periodic replay needs the same id for the same state every
period). In\--process only.

@<Globals@>=
int hashagg;
static unsigned char* ha_keys; static long ha_n, ha_kcap;   /* distinct keys, insertion order */
static long* ha_slot; static long ha_nslot;                 /* open\--addressing slots (id+1; 0=empty) */
static u64* ha_ea; static long ha_ea_n, ha_ea_cap;          /* edge accumulator: (dst<<32)|src */
static unsigned char* ha_sb; static int ha_sl;              /* sort context for the id permutation */

@ @<Subroutines@>=
static int vput(unsigned char*p,u64 x); static u64 vget(unsigned char**pp); static int vlen(u64 x);   /* LEB128 varints (defined with the compressed tables); used here to encode/decode edges */
static u64 ha_hash(unsigned char*k,int len){ u64 h=1469598103934665603ULL; int i; for(i=0;i<len;i++){ h^=k[i]; h*=1099511628211ULL; } return h; }
static void ha_rehash(long ns){ ha_slot=xrealloc(ha_slot,ns*sizeof(long)); long i; for(i=0;i<ns;i++) ha_slot[i]=0; ha_nslot=ns;
  long mask=ns-1,id; for(id=0;id<ha_n;id++){ u64 h=ha_hash(ha_keys+id*(long)ha_sl,ha_sl); long p=h&mask; while(ha_slot[p]) p=(p+1)&mask; ha_slot[p]=id+1; } }
static long ha_get_or_add(unsigned char*key){ int len=g_keylen; ha_sl=len;
  if(ha_n*10 >= ha_nslot*7) ha_rehash(ha_nslot? ha_nslot*2 : 4096);
  long mask=ha_nslot-1; u64 h=ha_hash(key,len); long p=h&mask;
  while(ha_slot[p]){ long id=ha_slot[p]-1; if(memcmp(ha_keys+id*(long)len,key,len)==0) return id; p=(p+1)&mask; }
  long id=ha_n++; if((id+1)*(long)len>ha_kcap){ ha_kcap=ha_kcap*2+(long)len*4096; ha_keys=xrealloc(ha_keys,ha_kcap); }
  memcpy(ha_keys+id*(long)len,key,len); ha_slot[p]=id+1; return id; }
static inline void ha_ea_add(u64 e){ if(ha_ea_n>=ha_ea_cap){ ha_ea_cap=ha_ea_cap*2+65536; ha_ea=xrealloc(ha_ea,ha_ea_cap*sizeof(u64)); } ha_ea[ha_ea_n++]=e; }
static int ha_idxcmp(const void*A,const void*B){ long a=*(const long*)A,b=*(const long*)B; return memcmp(ha_sb+a*ha_sl, ha_sb+b*ha_sl, ha_sl); }
static int ha_u64cmp(const void*A,const void*B){ u64 a=*(const u64*)A,b=*(const u64*)B; return a<b?-1:a>b?1:0; }

@ @<Subroutines@>=
void emit(unsigned char*key,int len,long src){
  if(hashagg){ long id=ha_get_or_add(key); ha_ea_add(((u64)id<<32)|(u32)src); return; }
  int t=omp_get_thread_num();
  if(tnr[t]>=trcap[t]){trcap[t]=trcap[t]*2+65536;trc[t]=xrealloc(trc[t],trcap[t]*sizeof(Rec));}
  Rec*R=&trc[t][tnr[t]++]; R->src=(int)src; memcpy(R->key,key,len); }
@#
int rcmp_r(const void*A,const void*B,void*arg){ (void)arg; const Rec*a=A,*b=B;
  return memcmp(a->key,b->key,g_keylen); }
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
static int* rkl[NRNG]; static long rcnt[NRNG], rgcap[NRNG];
static Edge* re[NRNG]; static long rne[NRNG], recap[NRNG];
static long bnd[NTMAX][NRNG+1];
void sort_reduce(int s,int rec){
  int NT=omp_get_max_threads(), t, r; long i;
  g_keylen=qnew;
  if(g_keylen>KEYMAX){ fprintf(stderr,"KEYMAX=%d too small for qnew=%d (m=%d); raise KEYMAX and rebuild\n",KEYMAX,qnew,m); exit(4); }
#pragma omp parallel for schedule(dynamic,1)
  for(t=0;t<NT;t++) if(tnr[t]) qsort_r(trc[t],tnr[t],sizeof(Rec),rcmp_r,0);
  long tot=0; for(t=0;t<NT;t++) tot+=tnr[t];
  if(getenv("BUILDMEM")){
    fprintf(stderr,"  [mem s=%d] emit recs=%ld (%.2fGB Rec, inline key len %d) ; state ncur=%ld (%.2fGB)\n",
      s, tot, tot*(double)sizeof(Rec)/1073741824.0, g_keylen, ncur,
      ncur*(double)(sizeof(long)+sizeof(int))/1073741824.0); }
  int P=NT*8; if(P>NRNG) P=NRNG; if(P<1) P=1; if((long)P>tot && tot>0) P=(int)tot;
  int S=P*8; if((long)S>tot) S=(int)tot;
  static Samp* samp=0; static long sampcap=0; if(sampcap<S+1){ sampcap=S+1; samp=xrealloc(samp,sampcap*sizeof(Samp)); }
  long sc=0;
  for(t=0;t<NT && tot>0;t++){ long nt=tnr[t]; if(nt==0) continue; long take=(long)((double)S*nt/tot); if(take<1) take=1; long j;
    for(j=0;j<take && sc<S;j++){ long idx=(long)((double)j*nt/take); if(idx>=nt) idx=nt-1; Rec*R=&trc[t][idx]; samp[sc].k=R->key; samp[sc].len=g_keylen; sc++; } }
  S=(int)sc; qsort(samp,S,sizeof(Samp),sampcmp);
  static Samp spl[NRNG]; int np=0, p;
  for(p=1;p<P;p++){ long si=(long)((double)p*S/P); if(si>=S) si=S-1; if(S>0) spl[np++]=samp[si]; }
  for(t=0;t<NT;t++){ bnd[t][0]=0; bnd[t][P]=tnr[t];
    for(p=0;p<np;p++){ long lo=0,hi=tnr[t]; while(lo<hi){ long mid=(lo+hi)/2; Rec*R=&trc[t][mid];
      if(keycmp(R->key,g_keylen,spl[p].k,spl[p].len)<0) lo=mid+1; else hi=mid; } bnd[t][p+1]=lo; } }
#pragma omp parallel for schedule(dynamic,1)
  for(r=0;r<P;r++){
    int hp[NTMAX]; long hpos[NTMAX]; int hn=0, tt;
    for(tt=0;tt<NT;tt++) hpos[tt]=bnd[tt][r];
#define RLESS(a,b) ({ Rec*Ra=&trc[a][hpos[a]],*Rb=&trc[b][hpos[b]]; keycmp(Ra->key,g_keylen,Rb->key,g_keylen)<0; })
    for(tt=0;tt<NT;tt++) if(hpos[tt]<bnd[tt][r+1]){ int c=hn++; hp[c]=tt;
      while(c>0){ int pp=(c-1)/2; if(RLESS(hp[c],hp[pp])){ int x=hp[c];hp[c]=hp[pp];hp[pp]=x; c=pp; } else break; } }
    rkuse[r]=0; rcnt[r]=0; rne[r]=0; int curlen=-1; long curpos=0;
    while(hn>0){ int mt=hp[0]; Rec*R=&trc[mt][hpos[mt]]; unsigned char*rkey=R->key; int rlen=g_keylen;
      int same=(curlen==rlen && memcmp(rk[r]+curpos,rkey,rlen)==0);
      if(!same){ if(curlen>=0){ if(rcnt[r]>=rgcap[r]){ rgcap[r]=rgcap[r]*2+1024; rkl[r]=xrealloc(rkl[r],rgcap[r]*sizeof(int)); } rkl[r][rcnt[r]]=curlen; rcnt[r]++; }
        if(rkuse[r]+rlen>rkcap[r]){ rkcap[r]=rkcap[r]*2+rlen+4096; rk[r]=xrealloc(rk[r],rkcap[r]); }
        curpos=rkuse[r]; memcpy(rk[r]+rkuse[r],rkey,rlen); rkuse[r]+=rlen; curlen=rlen; }
      if(rec){ if(rne[r]>=recap[r]){ recap[r]=recap[r]*2+1024; re[r]=xrealloc(re[r],recap[r]*sizeof(Edge)); } re[r][rne[r]].src=(int)R->src; re[r][rne[r]].dst=(int)rcnt[r]; re[r][rne[r]].c=1; rne[r]++; }
      hpos[mt]++;
      if(hpos[mt]>=bnd[mt][r+1]) hp[0]=hp[--hn];
      { int c=0; while(1){ int l=2*c+1,r2=2*c+2,sm=c; if(l<hn&&RLESS(hp[l],hp[sm]))sm=l; if(r2<hn&&RLESS(hp[r2],hp[sm]))sm=r2; if(sm==c)break; int x=hp[c];hp[c]=hp[sm];hp[sm]=x; c=sm; } }
    }
    if(curlen>=0){ if(rcnt[r]>=rgcap[r]){ rgcap[r]=rgcap[r]*2+1024; rkl[r]=xrealloc(rkl[r],rgcap[r]*sizeof(int)); } rkl[r][rcnt[r]]=curlen; rcnt[r]++; }
  }
  static long base[NRNG+1], kbase[NRNG+1], ebase[NRNG+1];
  base[0]=kbase[0]=ebase[0]=0;
  for(r=0;r<P;r++){ base[r+1]=base[r]+rcnt[r]; kbase[r+1]=kbase[r]+rkuse[r]; ebase[r+1]=ebase[r]+rne[r]; }
  ncur=base[P];
  curoff=xrealloc(curoff,(ncur+1)*sizeof(long)); curkl=xrealloc(curkl,(ncur+1)*sizeof(int)); curkp=xrealloc(curkp,kbase[P]+1);
#pragma omp parallel for schedule(dynamic,1)
  for(r=0;r<P;r++){ long off=kbase[r], j; memcpy(curkp+kbase[r], rk[r], rkuse[r]);
    for(j=0;j<rcnt[r];j++){ curkl[base[r]+j]=rkl[r][j]; curoff[base[r]+j]=off; off+=rkl[r][j]; } }
  for(t=0;t<NT;t++){ tnr[t]=0; }
  if(rec){ int r2; static long rne2[NRNG];
    /* Coalesce each range's edges IN PLACE and globalize dst -- parallel, no
       global consolidation buffer. All edges to a given dst live in one range
       (its key does), so per-range coalescing equals global coalescing; and the
       ranges hold disjoint, dst-ordered blocks, so concatenating them is already
       globally (dst,src)-sorted -- byte-identical to a global sort+coalesce, at
       half the peak edge memory (no |eb| copy) and cheaper (P small sorts). */
#pragma omp parallel for schedule(dynamic,1)
    for(r2=0;r2<P;r2++){ long pp,o2=0;
      if(rne[r2]>1) qsort(re[r2],rne[r2],sizeof(Edge),ecmp);
      for(pp=0;pp<rne[r2];){ long q2=pp+1; u64 cc=1;
        while(q2<rne[r2] && re[r2][q2].src==re[r2][pp].src && re[r2][q2].dst==re[r2][pp].dst){cc++;q2++;}
        re[r2][o2]=re[r2][pp]; re[r2][o2].c=cc; re[r2][o2].dst += base[r2]; o2++; pp=q2; }
      rne2[r2]=o2; }
    long o=0; for(r2=0;r2<P;r2++) o+=rne2[r2];
    nstate[reclev]=in_ncur_g; nstate[reclev+1]=ncur; nedge[reclev]=o; ncomp[reclev]=reccomp_n;
    if(stream_edges){  /* append this level to the open table file (prefix already written) */
      fwrite(&o,sizeof(long),1,rec_fp);
      for(r2=0;r2<P;r2++) if(rne2[r2]) fwrite(re[r2],sizeof(Edge),rne2[r2],rec_fp);
      fwrite(&reccomp_n,sizeof(long),1,rec_fp); fwrite(reccomp_buf,sizeof(Comp),reccomp_n,rec_fp);
      if(ferror(rec_fp)){ fprintf(stderr,"disk write error on %s (out of space?)\n",rec_path); exit(3); }
      edges[reclev]=0; comps[reclev]=0;
    } else {  /* in-process: compress the coalesced, dst-sorted edges (scheme-A varints)
                 straight into RAM -- ~4-6 B/edge instead of a raw 16 B Edge, so the
                 from-col-0 edge tables fit in memory. The ranges are disjoint
                 dst-ordered blocks, so we compress them in PARALLEL: pass 1 measures
                 each range's byte length (cross-range dst-delta included) to place it
                 in the final buffer, pass 2 encodes into that slot -- no serial
                 bottleneck and no memory spike. */
      double tc_=nowsec();
      static long rpd[NRNG], rolen[NRNG], roff[NRNG];
      for(r2=0;r2<P;r2++){ long pd=-1; int tt; for(tt=r2-1;tt>=0;tt--) if(rne2[tt]>0){ pd=re[tt][rne2[tt]-1].dst; break; } rpd[r2]=pd; }
#pragma omp parallel for schedule(dynamic,1)
      for(r2=0;r2<P;r2++){ long e2=0,len=0,prev_dst=rpd[r2];
        while(e2<rne2[r2]){ long dst=re[r2][e2].dst; len+=vlen((u64)(dst-prev_dst));
          long j=e2; while(j<rne2[r2] && re[r2][j].dst==dst) j++; len+=vlen((u64)(j-e2));
          long prev_src=0,k; for(k=e2;k<j;k++){ len+=vlen((u64)(re[r2][k].src-prev_src)); prev_src=re[r2][k].src; len+=vlen(re[r2][k].c); }
          prev_dst=dst; e2=j; } rolen[r2]=len; }
      long tot=0; for(r2=0;r2<P;r2++){ roff[r2]=tot; tot+=rolen[r2]; }
      unsigned char* ob=xmalloc(tot+16);
#pragma omp parallel for schedule(dynamic,1)
      for(r2=0;r2<P;r2++){ long e2=0,olen=roff[r2],prev_dst=rpd[r2];
        while(e2<rne2[r2]){ long dst=re[r2][e2].dst; olen+=vput(ob+olen,(u64)(dst-prev_dst));
          long j=e2; while(j<rne2[r2] && re[r2][j].dst==dst) j++; olen+=vput(ob+olen,(u64)(j-e2));
          long prev_src=0,k; for(k=e2;k<j;k++){ olen+=vput(ob+olen,(u64)(re[r2][k].src-prev_src)); prev_src=re[r2][k].src; olen+=vput(ob+olen,re[r2][k].c); }
          prev_dst=dst; e2=j; } }
      cedge[reclev]=ob; ced_len[reclev]=tot; edges[reclev]=0;
      /* the ranges are the parallel-decode split points: each starts at byte roff[r2]
         with running dst rpd[r2] and holds complete dst-groups (disjoint dst blocks),
         so a decoder can start there and its writes never collide with another range's. */
      { int sn=0,r3; csoff[reclev]=xrealloc(csoff[reclev],(P+1)*sizeof(long)); csdst[reclev]=xrealloc(csdst[reclev],(P+1)*sizeof(long));
        for(r3=0;r3<P;r3++) if(rne2[r3]>0){ csoff[reclev][sn]=roff[r3]; csdst[reclev][sn]=rpd[r3]; sn++; } csn[reclev]=sn; }
      t_compress+=nowsec()-tc_;
      comps[reclev]=xmalloc((reccomp_n+1)*sizeof(Comp)); memcpy(comps[reclev],reccomp_buf,reccomp_n*sizeof(Comp));
    }
    reclev++; }
}

@ |hashagg_finalize| plays the role of |sort_reduce| for the hash path: the
distinct keys are already deduplicated (in the hash), so we only sort them to a
canonical order, remap the accumulated edges' |dst| through that permutation, and
coalesce. In\--process only (no streaming yet).

@<Subroutines@>=
void hashagg_finalize(int rec){
  long N=ha_n, len=g_keylen, i;
  static long* idx=0; static long idxcap=0; if(idxcap<N+1){ idxcap=N+1; idx=xrealloc(idx,idxcap*sizeof(long)); }
  for(i=0;i<N;i++) idx[i]=i;
  ha_sb=ha_keys; ha_sl=(int)len; qsort(idx,N,sizeof(long),ha_idxcmp);
  static int* perm=0; static long permcap=0; if(permcap<N+1){ permcap=N+1; perm=xrealloc(perm,permcap*sizeof(int)); }
  for(i=0;i<N;i++) perm[idx[i]]=(int)i;
  ncur=N; curoff=xrealloc(curoff,(N+1)*sizeof(long)); curkl=xrealloc(curkl,(N+1)*sizeof(int)); curkp=xrealloc(curkp,N*len+1);
  for(i=0;i<N;i++){ memcpy(curkp+i*len, ha_keys+idx[i]*len, len); curoff[i]=i*len; curkl[i]=(int)len; }
  long o=0;
  if(rec){
    for(i=0;i<ha_ea_n;i++){ u64 e=ha_ea[i]; long od=(long)(e>>32); u32 sr=(u32)e; ha_ea[i]=((u64)(u32)perm[od]<<32)|sr; }
    qsort(ha_ea,ha_ea_n,sizeof(u64),ha_u64cmp);
    /* Encode scheme-A varints DIRECTLY from the sorted accumulator -- coalesce on
       the fly, no raw 16B-Edge intermediate (that was the finalize memory spike).
       Buffer grows from ~8B/edge; the compressed edges then live in RAM at
       ~6B/edge instead of 16B. */
    long cap=ha_ea_n*8+4096; unsigned char* ob=xmalloc(cap); long olen=0, prev_dst=-1, pp=0;
    long SPLIT=1L<<21, since=SPLIT, sn=0, scap=ha_ea_n/SPLIT+4;   /* split points every ~SPLIT edges for parallel decode */
    csoff[reclev]=xrealloc(csoff[reclev],scap*sizeof(long)); csdst[reclev]=xrealloc(csdst[reclev],scap*sizeof(long));
    while(pp<ha_ea_n){
      if(olen+64>cap){ cap=cap*2; ob=xrealloc(ob,cap); }
      if(since>=SPLIT){ csoff[reclev][sn]=olen; csdst[reclev][sn]=prev_dst; sn++; since=0; }
      u32 dd=(u32)(ha_ea[pp]>>32); long j=pp; while(j<ha_ea_n && (u32)(ha_ea[j]>>32)==dd) j++; since+=j-pp;
      long gsize=0,k=pp; while(k<j){ u32 s=(u32)ha_ea[k]; long q=k; while(q<j&&(u32)ha_ea[q]==s)q++; gsize++; k=q; }
      olen+=vput(ob+olen,(u64)((long)dd-prev_dst)); olen+=vput(ob+olen,(u64)gsize);
      long prev_src=0; k=pp;
      while(k<j){ if(olen+32>cap){ cap=cap*2; ob=xrealloc(ob,cap); } u32 s=(u32)ha_ea[k]; long q=k; u64 c=0; while(q<j&&(u32)ha_ea[q]==s){c++;q++;}
        olen+=vput(ob+olen,(u64)((long)s-prev_src)); prev_src=s; olen+=vput(ob+olen,c); o++; k=q; }
      prev_dst=dd; pp=j;
    }
    cedge[reclev]=ob; ced_len[reclev]=olen; edges[reclev]=0; nedge[reclev]=o; csn[reclev]=(int)sn;
    nstate[reclev]=in_ncur_g; nstate[reclev+1]=ncur; ncomp[reclev]=reccomp_n;
    comps[reclev]=xmalloc((reccomp_n+1)*sizeof(Comp)); memcpy(comps[reclev],reccomp_buf,reccomp_n*sizeof(Comp));
    reclev++;
  }
  if(getenv("BUILDMEM")) fprintf(stderr,"  [hamem] ea=%ld (%.2fGB 8B) -> comp %.2fGB (%.1fB/edge) + keys %.2fGB ; state=%ld\n",
    ha_ea_n, ha_ea_n*8.0/1073741824.0, rec?ced_len[reclev-1]/1073741824.0:0.0, (rec&&o)?(double)ced_len[reclev-1]/o:0.0, ha_n*(double)len/1073741824.0, ncur);
  ha_n=0; ha_ea_n=0; { long z; for(z=0;z<ha_nslot;z++) ha_slot[z]=0; }
}

@* The transfer step.
The bucket holds states as (key,weight) over the previous frontier. To place
cell~$s$: build the base mate table |bmate| over the new frontier (survivors
carried, newcomers bare) with $s$ in a temporary slot, then bring $s$ to
degree~2 by adding its $2-\deg(s)$ edges to future neighbours or the apex. When
|rec| is set we also record the integer edge table (src index $\to$ dst index)
and the completion contributions for this level.

@<Globals@>=
unsigned char*curkp; long*curoff; int*curkl; long ncur;
u128 cnt[1<<16];
u64 MODP=0;   /* 0 = exact u64; set to a prime for one CRT residue */
static inline u128 red(u128 x){ return MODP? x%MODP : x; }
int qnew, posS, apexnew, STEMP; int bmate[MAXF], o2n[MAXF];
#pragma omp threadprivate(bmate,STEMP)
@#
/* Reflection quotient. |sym_at[p]| is 1 when, after processing the in\--column
   position~$p$, the processed row set $\{|rowat|[0..p]\}$ is closed under the
   board reflection $r\to m-1-r$ (so the full processed region is symmetric and
   the reflection maps states to states). At such a cell we canonicalise each
   emitted key to $\min(k,\hat k)$ under the frontier\--position permutation
   |sym_sigma|; equal (also column\--shifted) reflected patterns then merge in
   |sort_reduce|. The DP runs in orbit\--sum coordinates: each surviving state
   carries the sum of its orbit's weights, expansion stays one\--rep\--per\--orbit,
   and completions credit that summed weight -- so counts are unchanged while the
   symmetric levels shrink by up to $2\times$. Valid only at symmetric cells;
   between them the sweep is ordinary. */
int sym_fold_on, sym_at[64], sym_bnd, sym_sigma[MAXF];
int recording, reclev;
typedef struct{ int src,dst; u64 c; } Edge;
typedef struct{ int src,delta; u64 mult; } Comp;
Edge* edges[MAXLEV]; long nedge[MAXLEV];
Comp* comps[MAXLEV]; long ncomp[MAXLEV];
long nstate[MAXLEV+1]; int Plevs;
int wfree, Ltot_g;   /* weight-free replay: record edge tables from col 0, iterate SpMV
                        from the trivial start [1] (regenerates all weights) instead of
                        capturing a seed at c0. Ltot_g = total recorded levels (cols 0..c0+period). */
int no_replay;       /* build only (leave the edge tables in RAM; the caller replays) */
int exact_capture;   /* the replay fills cnt2 but prints nothing (residues read by the exact CRT driver) */
static double nowsec(void);
double t_expand, t_reduce, t_compress;   /* build-phase timers (seconds) */
long del_total, del_total2;   /* DELCHECK: states that vanished vs same-phase prev col (period 1 / period 2) */
long memo_total; unsigned char* memo_buf[64]; long memo_bn[64], memo_cap[64]; int memo_kl[64];   /* MEMOSTAT: re-expansion factor */
/* Aggressive recording: record only cols 0..aggr_R (the left-boundary transient)
   plus the periodic block, SKIPPING the middle cols R+1..c0. Every column's state
   set is a subset of the stable set (verified), and the transfer is bulk from a
   couple columns in, so the periodic block reconstructs the middle columns. At the
   junction we remap the col-R boundary vector into the stable set's canonical order
   (|junc_remap[i]| = stable index of col-R state i, matched by key). -1 = off. */
int aggr_R=-1, junc_P0;   /* junc_P0 = #prefix levels = (aggr_R+1)*m */
unsigned char* junc_kr; long junc_nkr; int junc_klen;   /* col-R boundary keys */
int* junc_remap; long junc_nstab;                       /* remap[i] -> stable index; #stable states */
int period_g, c0_g, recend_col_g;   /* saved for dumping the extracted tables */
int direct_valid_col_g;             /* direct |cnt| is exact for columns $\le$ this */
int stop_after_record=0; int direct_cov_g=0;
int stream_edges=0; FILE* rec_fp=0; char rec_path[4096]; long nstate_off=0;  /* out-of-core recording */
int stream_spmv=0; FILE* tbl_fp=0; long lev_edge_off[MAXLEV];  /* out-of-core SpMV (edges streamed from disk) */
static u64 colfp_g[4096]; const char* ckpt_path=0; int ckpt_every=4;   /* checkpoint every K columns */
static unsigned char* ksbuf[3]; static long kscnt[3], kslen[3];   /* rolling exact boundary key sets (last 3 cols) for period detection */
u128* seedv; long seedn;
Comp* reccomp_buf; long reccomp_n, reccomp_cap;
int ecmp(const void*A,const void*B){ const Edge*a=A,*b=B; if(a->dst!=b->dst) return a->dst-b->dst; return a->src-b->src; }

@ @<Subroutines@>=
void build_bmate(int*omate,int qold){ int i;
  for(i=0;i<=qnew;i++) bmate[i]=-2; STEMP=qnew;
  for(i=0;i<qold;i++){ int dst=(i==posS)?STEMP:o2n[i]; if(dst<0) continue;
    if(omate[i]==-1) bmate[dst]=-1;
    else if(omate[i]>=0){ int op=omate[i]; int pdst=(op==posS)?STEMP:o2n[op]; bmate[dst]=pdst>=0?pdst:-2; } } }

@ |canon_nk| replaces a just\--built key by the lexicographically smaller of it
and its reflection (identity unless |sym_bnd| is set for this cell). It reads
only |sym_bnd| and |sym_sigma|, fixed before the parallel expansion, so it is
safe to call from each thread.

@<Subroutines@>=
static inline void canon_nk(unsigned char*nk,int nl){
  if(!sym_bnd) return;
  unsigned char rk[MAXF]; int i;
  for(i=0;i<nl;i++){ int b=nk[i], d=sym_sigma[i];
    rk[d] = b==0?0 : b==255?255 : (unsigned char)(1+sym_sigma[b-1]); }
  if(memcmp(rk,nk,nl)<0) memcpy(nk,rk,nl); }

@ @<Subroutines@>=
void run_step(int s,int rec){
  int i; long si; long in_ncur=ncur; int qold=frontier_before(s); posS=ifrb[s];
  static int frold[MAXF]; for(i=0;i<qold;i++) frold[i]=fr[i];
  qnew=frontier_before(s+1); static int ifrnew[64*64+2], frnew[MAXF];
  for(i=0;i<qnew;i++) frnew[i]=fr[i]; for(i=0;i<=V;i++) ifrnew[i]=-1; for(i=0;i<qnew;i++) ifrnew[frnew[i]]=i;
  apexnew=ifrnew[V];
  sym_bnd = sym_fold_on && sym_at[s%m];
  if(sym_bnd){ for(i=0;i<qnew;i++){ int v=frnew[i], rv;
      if(v==V) rv=V; else { int col=v/m, pp=v%m; rv=col*m + posof[m-1-rowat[pp]]; }
      sym_sigma[i]=ifrnew[rv]; } }
  int nbr[16],rr=0; nbr[rr++]=apexnew;
  for(i=0;i<ND[s];i++){ int w=NB[s][i]; if(w>s && ifrnew[w]>=0) nbr[rr++]=ifrnew[w]; }
  for(i=0;i<qold;i++) o2n[i]=(frold[i]==s)?-1:ifrnew[frold[i]];
  in_ncur_g=ncur; reccomp_n=0; g_keylen=qnew;
  { double te_=nowsec(); @<Expand every state at cell |s|@>@; t_expand+=nowsec()-te_; }
  { double tr_=nowsec(); @<Sort, reduce, and (if recording) capture the edge table@>@; t_reduce+=nowsec()-tr_; }
}

@ Each state expands independently, so the loop runs in parallel (except while
recording, when it stays serial to keep the completion log ordered). The
transition scratch (|mate|, |bmate|, |cycle|, |STEMP|) is |threadprivate|; each
thread emits into its own pool; |cnt| is updated in a critical section.

@<Expand every state at cell |s|@>=
#pragma omp parallel for schedule(dynamic,32) if(!hashagg)
for(si=0;si<ncur;si++){
  int i,a,b,nl,deg,need,mp; unsigned char*ok=curkp+curoff[si]; int okl=curkl[si];
  int omate[MAXF]; unsigned char nk[MAXF];
  for(i=0;i<okl;i++){int cc=ok[i]; omate[i]= cc==0?-2 : cc==255?-1 : cc-1;}
  deg=omate[posS]==-2?0:omate[posS]==-1?2:1; need=2-deg;
  if(need==0){ build_bmate(omate,qold); for(i=0;i<qnew;i++) mate[i]=bmate[i]; nl=keyof(qnew,nk); canon_nk(nk,nl); emit(nk,nl,si); }
  else if(need==1){ for(a=0;a<rr;a++){ build_bmate(omate,qold); for(i=0;i<=qnew;i++) mate[i]=bmate[i];
    if(add_derived(STEMP,nbr[a])){ nl=keyof(qnew,nk); canon_nk(nk,nl); emit(nk,nl,si); }
    else if(cycle){ mp=completion_mp(qnew,apexnew,frnew,s); if(mp) @<Credit completion@>; } } }
  else { for(a=0;a<rr;a++)for(b=a+1;b<rr;b++){ build_bmate(omate,qold); for(i=0;i<=qnew;i++) mate[i]=bmate[i];
    if(!add_derived(STEMP,nbr[a])){ if(cycle){ mp=completion_mp(qnew,apexnew,frnew,s); if(mp) @<Credit completion@>; } continue; }
    if(add_derived(STEMP,nbr[b])){ nl=keyof(qnew,nk); canon_nk(nk,nl); emit(nk,nl,si); }
    else if(cycle){ mp=completion_mp(qnew,apexnew,frnew,s); if(mp) @<Credit completion@>; } } }
}

@ @<Credit completion@>=
{
#pragma omp critical
  { if(rec){ if(reccomp_n>=reccomp_cap){ reccomp_cap=reccomp_cap*2+(1<<20); reccomp_buf=xrealloc(reccomp_buf,reccomp_cap*sizeof(Comp)); }
      reccomp_buf[reccomp_n].src=si; reccomp_buf[reccomp_n].delta=mp-(s+1); reccomp_buf[reccomp_n].mult=1; reccomp_n++; } }
}

@ The sort and reduce run across the per-thread pools in parallel (see
|sort_reduce|); |in_ncur_g| records the input level size for the edge table.

@<Sort, reduce, and (if recording) capture the edge table@>=
if(hashagg) hashagg_finalize(rec); else sort_reduce(s,rec);

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
  fwrite(cnt,sizeof(u128),1<<16,f); fwrite(colfp_g,sizeof(u64),4096,f);
  fwrite(&ncur,sizeof(long),1,f); fwrite(&kuse,sizeof(long),1,f);
  fwrite(curoff,sizeof(long),ncur,f); fwrite(curkl,sizeof(int),ncur,f);
  /* weight-free: no state weights to checkpoint */
  { long kb=0,i; for(i=0;i<ncur;i++) kb+=curkl[i]; fwrite(&kb,sizeof(long),1,f); fwrite(curkp,1,kb,f); }
  fclose(f); }
int load_ckpt(const char*path){ FILE*f=fopen(path,"rb"); if(!f) return -1; int nexts,mm;
  if(fread(&mm,sizeof(int),1,f)!=1) return -1; m=mm; fread(&nexts,sizeof(int),1,f);
  fread(cnt,sizeof(u128),1<<16,f); fread(colfp_g,sizeof(u64),4096,f);
  fread(&ncur,sizeof(long),1,f); fread(&kuse,sizeof(long),1,f);
  curoff=xrealloc(curoff,(ncur+1)*sizeof(long)); curkl=xrealloc(curkl,(ncur+1)*sizeof(int));
  fread(curoff,sizeof(long),ncur,f); fread(curkl,sizeof(int),ncur,f);
  { long kb; fread(&kb,sizeof(long),1,f); curkp=xrealloc(curkp,kb+1); fread(curkp,1,kb,f); }
  fclose(f); fprintf(stderr,"resumed build from %s at cell %d (col %d)\n",path,nexts,nexts/mm); return nexts; }

@ @* Dumping and reloading the extracted tables.
The expensive part is the build that discovers and records the periodic transfer.
Once recorded, the tables (edge tables, completions, seed vector, level sizes,
period, and the directly\--counted columns up to |recend_col_g|) are small and
self\--contained, so we dump them to disk. A later run reloads them and does only
the cheap SpMV --- and a build that dies can be re\--run without touching the SpMV.

@<Subroutines@>=
void dump_tables(const char*path){ int L;
  if(stream_edges && rec_fp){
    /* The prefix (hdr, seed, nstate placeholder) and every level are already on
       disk (written during recording). Append cnt, then patch the two fields not
       known until now: hdr[5]=direct_valid_col_g and the real nstate[]. No copy,
       no second file -- so the table needs disk for its own size only. */
    int nc=(direct_valid_col_g+1)*m; fwrite(&nc,sizeof(int),1,rec_fp); fwrite(cnt,sizeof(u128),nc,rec_fp);
    fflush(rec_fp);
    fseek(rec_fp,5L*sizeof(int),SEEK_SET); fwrite(&direct_valid_col_g,sizeof(int),1,rec_fp);
    fseek(rec_fp,nstate_off,SEEK_SET); fwrite(nstate,sizeof(long),Plevs+1,rec_fp);
    if(ferror(rec_fp)){ fprintf(stderr,"disk write error on %s (out of space?)\n",path); exit(3); }
    fclose(rec_fp); rec_fp=0;
    fprintf(stderr,"dumped tables to %s (period %d, c0 %d, %d levels)\n",path,period_g,c0_g,Plevs); return; }
  FILE*f=fopen(path,"wb");   /* non-streaming fallback (edges held in RAM) */
  int hdr[6]={m,period_g,c0_g,Plevs,recend_col_g,direct_valid_col_g}; fwrite(hdr,sizeof(int),6,f);
  fwrite(&seedn,sizeof(long),1,f); fwrite(seedv,sizeof(u128),seedn,f);
  fwrite(nstate,sizeof(long),Plevs+1,f);
  for(L=0;L<Plevs;L++){ fwrite(&nedge[L],sizeof(long),1,f); fwrite(edges[L],sizeof(Edge),nedge[L],f);
    fwrite(&ncomp[L],sizeof(long),1,f); fwrite(comps[L],sizeof(Comp),ncomp[L],f); }
  int nc=(direct_valid_col_g+1)*m; fwrite(&nc,sizeof(int),1,f); fwrite(cnt,sizeof(u128),nc,f);
  fclose(f); fprintf(stderr,"dumped tables to %s (period %d, c0 %d, %d levels)\n",path,period_g,c0_g,Plevs); }
int load_tables(const char*path){ FILE*f=fopen(path,"rb"); if(!f) return 0; int L,hdr[6];
  if(fread(hdr,sizeof(int),6,f)!=6) return 0; m=hdr[0]; period_g=hdr[1]; c0_g=hdr[2]; Plevs=hdr[3]; recend_col_g=hdr[4]; direct_valid_col_g=hdr[5];
  fread(&seedn,sizeof(long),1,f); seedv=xmalloc(seedn*sizeof(u128)); fread(seedv,sizeof(u128),seedn,f);
  fread(nstate,sizeof(long),Plevs+1,f);
  /* Stream edges from disk (tables can be >RAM, e.g. 100GB for m=7): keep the
     small comps[] in RAM, remember each level's edge offset, skip the edges. */
  for(L=0;L<Plevs;L++){ fread(&nedge[L],sizeof(long),1,f);
    lev_edge_off[L]=ftello(f); fseeko(f,(off_t)nedge[L]*sizeof(Edge),SEEK_CUR);
    fread(&ncomp[L],sizeof(long),1,f); comps[L]=xmalloc((ncomp[L]+1)*sizeof(Comp)); fread(comps[L],sizeof(Comp),ncomp[L],f);
    edges[L]=0; }
  int nc; fread(&nc,sizeof(int),1,f); memset(cnt,0,sizeof(cnt)); fread(cnt,sizeof(u128),nc,f);
  tbl_fp=f; stream_spmv=1; return 1; }   /* keep file open for edge streaming */

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
  curkp=xmalloc(1<<20); curoff=xmalloc(sizeof(long)); curkl=xmalloc(sizeof(int));
  memcpy(curkp,key,kl); curoff[0]=0; curkl[0]=kl; ncur=1; }
@#
u64 cnt2[1<<16];
int run_periodic(int mm,int Wb,int Next){
  int i,s,c; m=mm; n=Wb; build_board(); memset(cnt2,0,sizeof(cnt2));
  sym_fold_on = getenv("SYMFOLD")?1:0;
  if(sym_fold_on) fprintf(stderr,"reflection fold ON (rowat outside-in? %s)\n", getenv("OUTIN")?"yes":"no");
  wfree = 1; Ltot_g=0;   /* weight-free replay is the only build mode now */
  t_expand=t_reduce=t_compress=0; del_total=del_total2=0;
  if(getenv("MEMOSTAT")){ memo_total=0; int ss; for(ss=0;ss<64;ss++) memo_bn[ss]=0; }
  hashagg = getenv("HASHAGG")?1:0;
  aggr_R = getenv("AGGR")? atoi(getenv("AGGR")) : -1;   /* aggressive: record only cols 0..AGGR + periodic block */
  if(hashagg && stream_edges){ fprintf(stderr,"HASHAGG is in-process only for now\n"); exit(4); }
  { int L; for(L=0;L<MAXLEV;L++){ if(cedge[L]){free(cedge[L]);cedge[L]=0;} edges[L]=0; ced_len[L]=0; nedge[L]=0;
      if(csoff[L]){free(csoff[L]);csoff[L]=0;} if(csdst[L]){free(csdst[L]);csdst[L]=0;} csn[L]=0; } }  /* fresh run: clear per-level edge stores + split points */
  recording=0; reclev=0; c0_g=-1; recend_col_g=-1; direct_valid_col_g=-1; Plevs=0; seedn=0;
  if(resume_from<=0){ memset(cnt,0,sizeof(cnt)); seed_bucket(); }  /* fresh; else state is loaded */
  int c0=-1,period=0,recstart=-1,recend=-1,seam_break=-1;
  int cand_p=0,cand_c=-1;   /* pending (unconfirmed) period candidate */
  if(wfree){ recording=1; recstart=0; recend=V; }  /* record every level from col 0; recend widened until the period is confirmed */
  for(s=(resume_from>0?resume_from:0);s<V;s++){
    if(recording && s==recstart){ seedn=ncur; nstate[0]=ncur; reclev=0;  /* weight-free: no seed vector; start replay from [1] */
      if(getenv("SYMDUMP")){  /* dump frontier layout + all state keys, for the reflection-symmetry study */
        int q=frontier_before(recstart); FILE*sp=fopen(getenv("SYMDUMP"),"wb");
        fwrite(&m,sizeof(int),1,sp); fwrite(&V,sizeof(int),1,sp); fwrite(&q,sizeof(int),1,sp); fwrite(fr,sizeof(int),q,sp);
        fwrite(&ncur,sizeof(long),1,sp); { long i2; for(i2=0;i2<ncur;i2++){ fwrite(&curkl[i2],sizeof(int),1,sp); fwrite(curkp+curoff[i2],1,curkl[i2],sp); } }
        fclose(sp); fprintf(stderr,"symdump: m=%d V=%d q=%d ncur=%ld -> %s\n",m,V,q,ncur,getenv("SYMDUMP")); }
      if(stream_edges){  /* write the table's fixed prefix now; levels stream in after; patch at dump */
        rec_fp=fopen(rec_path,"wb"); if(!rec_fp){ fprintf(stderr,"cannot open %s\n",rec_path); exit(3); }
        int hdr[6]={m,period_g,c0_g,Plevs,recend_col_g,0};   /* direct_valid_col_g patched at dump */
        fwrite(hdr,sizeof(int),6,rec_fp); fwrite(&seedn,sizeof(long),1,rec_fp); fwrite(seedv,sizeof(u128),seedn,rec_fp);
        nstate_off=ftell(rec_fp); { long z[MAXLEV+1]; int L; for(L=0;L<=Plevs;L++) z[L]=0; fwrite(z,sizeof(long),Plevs+1,rec_fp); }
        if(ferror(rec_fp)){ fprintf(stderr,"disk write error on %s (out of space?)\n",rec_path); exit(3); } } }
    if(aggr_R>=0 && s==(aggr_R+1)*m){  /* capture the col-R boundary keys (col-R canonical order) */
      junc_nkr=ncur; junc_klen=curkl[0]; junc_kr=xrealloc(junc_kr,ncur*(long)junc_klen+1);
      long i2; for(i2=0;i2<ncur;i2++) memcpy(junc_kr+i2*junc_klen, curkp+curoff[i2], junc_klen); }
    if(aggr_R>=0 && c0>=0 && period>0){  /* stable set reached at the phase matching col-R: build remap col-R -> stable by key.
        The replay applies the periodic block from col R+1, which starts at periodic column pi=(R-c0) mod period,
        whose INPUT is the col-(c0+pi) boundary -- capture THAT phase so the phases line up. */
      int pi=((aggr_R-c0)%period+period)%period;
      if(s==(c0+pi+1)*m){
        junc_nstab=ncur; junc_remap=xrealloc(junc_remap,(junc_nkr+1)*sizeof(int)); long i2;
        for(i2=0;i2<junc_nkr;i2++){ long lo=0,hi=ncur; unsigned char*kk=junc_kr+i2*junc_klen;
          while(lo<hi){ long mid=(lo+hi)/2; int d=memcmp(curkp+curoff[mid],kk,junc_klen); if(d<0)lo=mid+1; else hi=mid; }
          junc_remap[i2]=(lo<ncur && memcmp(curkp+curoff[lo],kk,junc_klen)==0)? (int)lo : -1; }
        fprintf(stderr,"aggr junction: col-%d states %ld -> stable(phase %d) %ld (remap built)\n",aggr_R,junc_nkr,pi,junc_nstab); } }
    { int recq; if(aggr_R>=0) recq = recording && ( s < (aggr_R+1)*m || (c0>=0 && s>=(c0+1)*m && s<=recend) );
      else recq = recording && s>=recstart && s<=recend;
      run_step(s, recq); }
    if(s%m==m-1){ c=s/m;
      { int slot=c%3; long tot=0,i2; for(i2=0;i2<ncur;i2++) tot+=curkl[i2];   /* keep the exact boundary key set (sorted) for cheap count + exact period test */
        ksbuf[slot]=xrealloc(ksbuf[slot],tot+1); long off=0; for(i2=0;i2<ncur;i2++){ memcpy(ksbuf[slot]+off,curkp+curoff[i2],curkl[i2]); off+=curkl[i2]; }
        kscnt[slot]=ncur; kslen[slot]=tot; }
      @<Detect and confirm the periodic column@>@; }
    if(s%m==m-1 && getenv("DUMPCOLS")){ /* dump this column's sorted key set (subset study) */
      char pth[4200]; snprintf(pth,sizeof pth,"%s/col%d.bin",getenv("DUMPCOLS"),s/m);
      FILE*dp=fopen(pth,"wb"); long t2; int ql=frontier_before(s+1);
      fwrite(&ql,sizeof(int),1,dp); fwrite(&ncur,sizeof(long),1,dp);
      for(t2=0;t2<ncur;t2++) fwrite(curkp+curoff[t2],1,curkl[t2],dp); fclose(dp); }
    if(s%m==m-1 && getenv("DBG")) fprintf(stderr,"  col %d ncur=%ld rec=%d\n",s/m,ncur,recording);
    if(!recording && s%m==m-1 && (s/m)%ckpt_every==0) save_ckpt(s+1);
    if(recording && s==seam_break && stop_after_record) break;
    if(getenv("MEMOSTAT")){  /* A1 ceiling: total expansions vs DISTINCT (substep,state); ratio = re-expansion factor */
      int ss=s%m; long i2; memo_total += ncur;
      for(i2=0;i2<ncur;i2++){ int kl=curkl[i2]; if(memo_bn[ss]+kl>memo_cap[ss]){ memo_cap[ss]=memo_cap[ss]*2+kl+(1<<20); memo_buf[ss]=xrealloc(memo_buf[ss],memo_cap[ss]); }
        memcpy(memo_buf[ss]+memo_bn[ss], curkp+curoff[i2], kl); memo_bn[ss]+=kl; }
      memo_kl[ss]=curkl[0]; }
    if(getenv("DELCHECK")){  /* deletions vs SAME-PHASE previous column: s-m (period 1) and s-2m (period 2) */
      static unsigned char* cb[96]; static long cc[96]; static int ck[96]; int W=2*m+2, sl=s%W;
      long tot=0,i2; for(i2=0;i2<ncur;i2++) tot+=curkl[i2];
      int P; for(P=1;P<=2;P++){ int back=P*m; if(s<back) continue; int op=(s-back)%W; if(ck[op]!=curkl[0]) continue;
        int kl=ck[op]; long a,b=0,del=0;
        for(a=0;a<cc[op];a++){ unsigned char*pk=cb[op]+a*kl;
          while(b<ncur && memcmp(curkp+curoff[b],pk,kl)<0) b++;
          if(b<ncur && memcmp(curkp+curoff[b],pk,kl)==0) b++; else del++; }
        if(P==1) del_total+=del; else del_total2+=del;
        if(del && s%m==m-1) fprintf(stderr,"  col %d boundary: %ld states in col-%d(p=%d) not here\n",s/m,del,s/m-P,P); }
      cb[sl]=xrealloc(cb[sl],tot+1); { long off=0; for(i2=0;i2<ncur;i2++){memcpy(cb[sl]+off,curkp+curoff[i2],curkl[i2]); off+=curkl[i2];} } cc[sl]=ncur; ck[sl]=curkl[0]; }
    if(wfree && recording && recend>=0 && s>=recend) break;  /* wfree: nothing to sweep past the recorded block */
  }
  @<Finalize direct coverage and verify periodicity@>@;
  if(getenv("MEMOSTAT")){ long dist=0; int ss; for(ss=0;ss<m;ss++){ if(!memo_bn[ss])continue; int kl=memo_kl[ss]; long n=memo_bn[ss]/kl,i3;
      long* idx=xmalloc(n*sizeof(long)); for(i3=0;i3<n;i3++) idx[i3]=i3; ha_sb=memo_buf[ss]; ha_sl=kl; qsort(idx,n,sizeof(long),ha_idxcmp);
      long d=n?1:0; for(i3=1;i3<n;i3++) if(memcmp(memo_buf[ss]+idx[i3]*kl, memo_buf[ss]+idx[i3-1]*kl, kl)!=0) d++; dist+=d; free(idx); }
    fprintf(stderr,"MEMOSTAT: total expansions=%ld, distinct (substep,state)=%ld  => re-expansion factor %.2fx (A1 ceiling)\n",
      memo_total,dist, dist?(double)memo_total/dist:0.0); }
  if(getenv("DELCHECK")) fprintf(stderr,"DELCHECK: deletions vs prev col -- period1 (s-m)=%ld, period2 (s-2m)=%ld => %s\n",
    del_total, del_total2, (del_total==0||del_total2==0)?"MONOTONE at the true period (each state expandable once)":"deletions at both periods");
  direct_cov_g = stop_after_record ? direct_valid_col_g : n;
  if(!stream_edges && !no_replay) spmv_run(Next, 1);   /* streamed edges are freed; SpMV via a separate `run` */
  return c0;
}

@ We confirm the period {\bf exactly}, not by a fingerprint. At each boundary we
first compare the state {\it count\/} $|S_c|$ against $|S_{c-p}|$ (O(1)); only if
they match do we do the full byte compare of the sorted key sets ($S_c=S_{c-p}$).
An exact match is {\it definitive\/}: the transfer is deterministic, so identical
sets repeat forever -- no coincidence is possible (unlike a hash), and no
"persist one more period" wait is needed. We commit at the {\bf first} exact
repeat, which drops the conservative confirmation columns the old fingerprint
scheme swept. (The periodic block is still recorded ahead as cols $c_0+1\ldots
c_0+p$; since $S_{c_0}=S_{c_0-p}$ that is byte-identical to the preceding $p$
columns.)

@<Detect and confirm the periodic column@>=
if(c0<0){ int p;
  for(p=1;p<=2;p++){ if(c<p) continue; int cur=c%3, old=(c-p)%3;
    if(kscnt[cur]==kscnt[old] && kslen[cur]==kslen[old]
       && memcmp(ksbuf[cur],ksbuf[old],kslen[cur])==0 && c+p+3<=n){
      period=p; c0=c; recend=(c0+1+period)*m-1;
      Plevs=period*m; period_g=period; c0_g=c0; recend_col_g=c0+period;
      if(wfree){ seam_break=recend;
        if(aggr_R>=0){ junc_P0=(aggr_R+1)*m; Ltot_g=junc_P0+Plevs; }  /* prefix cols 0..R, then the periodic block */
        else Ltot_g=recend+1; }  /* recorded cols 0..c0+period */
      else { recstart=(c0+1)*m; recording=1; seam_break=recend+3*m; }
      fprintf(stderr,"stable col %d, period %d (exact match)%s\n",c0,period,
        wfree?(aggr_R>=0?" [wfree aggressive: col 0..R + periodic]":" [wfree: recorded from col 0]"):"");
      break; } }
}

@ The direct sweep counts every column it reaches, exact once that column's
completions have all closed (a lag of a couple knight\--moves). The SpMV cannot
reconstruct completions that close inside the recorded/seed columns, so its
{\it first\/} extrapolated column, |recend+1|, comes out short. We therefore
trust the direct count one column past the recorded region --- the |seam_break|
sweep ran far enough for its completions to close --- and hand over to the SpMV
only from |recend+2|, exactly where the cross\--check begins. If a period was
recorded we also assert the transfer is stationary; a mismatch means the strip
was too narrow and the recording is boundary\--contaminated.

@<Finalize direct coverage and verify periodicity@>=
{ int swept = (s>=V)? n : s/m;
  if(c0<0){ direct_valid_col_g = swept>0? swept-1 : 0;   /* no period: direct only */
    fprintf(stderr,"no confirmed period within strip W_b=%d: direct counts only "
      "(increase W_b to record the periodic transfer)\n",n); }
  else { direct_valid_col_g = (swept>=recend_col_g+3)? recend_col_g+1 : recend_col_g;
    /* stationarity: the periodic block must map its state set to an identical one.
       Non-wfree indexes the block at [0,Plevs]; wfree records from col 0, so its
       periodic block is the LAST Plevs levels: [Ltot-Plevs, Ltot]. */
    long sa = wfree? nstate[Ltot_g-Plevs] : nstate[0];
    long sb = wfree? nstate[Ltot_g]       : nstate[Plevs];
    if(sa!=sb){
      fprintf(stderr,"ERROR: recorded transfer not stationary (%ld vs %ld); strip "
        "W_b=%d too narrow -- recording columns are boundary-contaminated. "
        "Rebuild with a larger W_b.\n", sa,sb,n); exit(4); } }
}

@ @<Subroutines@>=
/* SpMV from the (built or loaded) tables out to column |Nto|. |crosscheck|
   compares against the direct |cnt| where both are valid. Uses |cnt2| scratch. */
int build_extract(int mm,int Wb){ return run_periodic(mm,Wb,Wb); }
void spmv_run(int Nto,int crosscheck){
  int i; memset(cnt2,0,sizeof(cnt2));
  @<Iterate the SpMV out to column $N$@>;
  { int cc;  /* weight-free: every column comes from the replay */
    if(!exact_capture) for(cc=1;cc<=Nto;cc++){ u64 val=cnt2[cc*m]; if(val) printf("open %dx%d = %llu\n",m,cc,val); }
    return; }
  { int good=1,cc;
    if(crosscheck) for(cc=c0_g+3;cc<=Nto && cc<=direct_cov_g;cc++) if((u64)cnt[cc*m]!=cnt2[cc*m]){ good=0;
      fprintf(stderr,"MISMATCH c=%d direct=%llu spmv=%llu\n",cc,(unsigned long long)(u64)cnt[cc*m],(unsigned long long)cnt2[cc*m]); }
    for(cc=1;cc<=Nto;cc++){ u64 val = cc<=direct_valid_col_g? (u64)red(cnt[cc*m]) : cnt2[cc*m]; if(val) printf("open %dx%d = %llu%s\n",m,cc,val, cc>direct_valid_col_g?"  (SpMV)":""); }
    if(crosscheck) fprintf(stderr,"%s\n", good?"SpMV matches direct":"SpMV MISMATCH"); }
}

@ @<Iterate the SpMV out to column $N$@>=
if(wfree && Ltot_g>0){
  /* Weight-free replay: start from the trivial vector at col 0 (a single state,
     weight 1), apply the recorded prefix (cols $0..c_0+p$), then repeat the
     periodic block. Regenerates every weight and completion; no seed needed. */
#define APPLY(LL,ACOL,SUB) do{ int L_=(LL),abscol_=(ACOL),substep_=(SUB); long e_; \
    for(e_=0;e_<ncomp[L_];e_++){ int idx_=abscol_*m+substep_+1+comps[L_][e_].delta; cnt2[idx_]=red(cnt2[idx_]+(u128)v[comps[L_][e_].src]*comps[L_][e_].mult); } \
    u64* vn_=xcalloc(nstate[L_+1],sizeof(u64)); \
    if(cedge[L_]){ unsigned char* p_=cedge[L_], *end_=p_+ced_len[L_]; long dst_=-1; \
      while(p_<end_){ dst_+=(long)vget(&p_); long gs_=(long)vget(&p_),g_,src_=0; \
        for(g_=0;g_<gs_;g_++){ src_+=(long)vget(&p_); u64 c_=vget(&p_); vn_[dst_]=red(vn_[dst_]+(u128)v[src_]*c_); } } } \
    else { long ne_=nedge[L_],e2_; for(e2_=0;e2_<ne_;e2_++){ int d_=edges[L_][e2_].dst; vn_[d_]=red(vn_[d_]+(u128)v[edges[L_][e2_].src]*edges[L_][e2_].c); } } \
    free(v); v=vn_; }while(0)
  u64* v=xcalloc(nstate[0],sizeof(u64)); v[0]=1;
  int L, X, ss, pfx=(aggr_R>=0)?junc_P0:Ltot_g, perbase=(aggr_R>=0)?junc_P0:(Ltot_g-Plevs), pref=c0_g+1;
  for(L=0;L<pfx;L++) APPLY(L, L/m, L%m);              /* prefix: cols 0..R (or 0..recend_col) */
  int basecol;
  if(aggr_R>=0){ u64* vs=xcalloc(junc_nstab,sizeof(u64)); long i2;   /* remap col-R -> stable order */
    for(i2=0;i2<junc_nkr;i2++){ int j=junc_remap[i2]; if(j>=0) vs[j]=v[i2]; }
    free(v); v=vs; basecol=aggr_R+1; }
  else basecol=recend_col_g+1;
  for(X=basecol; X<=Nto; X++){ int pc=((X-pref)%period_g+period_g)%period_g;   /* phase-aware: col X uses periodic column pc */
    for(ss=0;ss<m;ss++) APPLY(perbase+pc*m+ss, X, ss); }
  free(v);
#undef APPLY
} else
if(c0_g>=0 && Plevs>0){ u64* v=xmalloc(nstate[0]*sizeof(u64)); { long i2; for(i2=0;i2<nstate[0];i2++) v[i2]=(u64)(MODP?seedv[i2]%MODP:seedv[i2]); }
  int basecol=c0_g+1;
  while(basecol<=Nto){ int L;
    for(L=0;L<Plevs;L++){ int abscol=basecol+L/m, substep=L%m; long e;
      for(e=0;e<ncomp[L];e++){ int idx=abscol*m+substep+1+comps[L][e].delta; cnt2[idx]=red(cnt2[idx]+v[comps[L][e].src]*comps[L][e].mult); }
      u64* vn=xcalloc(nstate[L+1],sizeof(u64));
      { long ne=nedge[L];
        if(stream_spmv){  /* out-of-core: stream this level's edges from disk in chunks */
          fseeko(tbl_fp,lev_edge_off[L],SEEK_SET);
          static Edge* buf=0; static long bufcap=0; long CH=1L<<20;
          if(bufcap<CH){ bufcap=CH; buf=xrealloc(buf,bufcap*sizeof(Edge)); }
          long done=0;
          while(done<ne){ long chunk=ne-done; if(chunk>CH) chunk=CH;
            if((long)fread(buf,sizeof(Edge),chunk,tbl_fp)!=chunk){ fprintf(stderr,"table edge read error\n"); exit(3); }
            long e2; for(e2=0;e2<chunk;e2++){ int d=buf[e2].dst; vn[d]=red(vn[d]+v[buf[e2].src]*buf[e2].c); }
            done+=chunk; }
        } else {          /* in-RAM: dst-partitioned parallel apply */
          int NT=omp_get_max_threads(),tt; static long spl[NTMAX+1];
          spl[0]=0; spl[NT]=ne;
          for(tt=1;tt<NT;tt++){ long g=(long)((double)tt*ne/NT);
            while(g>0&&g<ne&&edges[L][g].dst==edges[L][g-1].dst) g++; spl[tt]=g; }
#pragma omp parallel for schedule(dynamic,1)
          for(tt=0;tt<NT;tt++){ long e2; for(e2=spl[tt];e2<spl[tt+1];e2++){ int d=edges[L][e2].dst;
            vn[d]=red(vn[d]+v[edges[L][e2].src]*edges[L][e2].c); } }
        } }
      free(v); v=vn; }
    basecol+=period_g; }
  free(v); }

@ @<Report and cross\--check@>=
{ int good=1,cc;
  for(cc=c0+3;cc<=Next && cc<=Wb;cc++) if(cnt[cc*m]!=cnt2[cc*m]){ good=0;
    fprintf(stderr,"MISMATCH c=%d direct=%llu spmv=%llu\n",cc,cnt[cc*m],cnt2[cc*m]); }
  for(cc=1;cc<=Next;cc++){ u64 val = cc<=Wb? cnt[cc*m] : cnt2[cc*m]; if(val) printf("open %dx%d = %llu%s\n",m,cc,val, cc>Wb?"  (SpMV)":""); }
  fprintf(stderr,"%s\n", good?"SpMV matches direct":"SpMV MISMATCH"); }

@ Block SpMV: compute several primes' residues in ONE pass over the edge table.
The per\--prime state vectors are interleaved (|v[i*B+b]|), so a cache line serves
several primes at once; the edges (sorted by |dst|) are summed per\--|dst| in a u64
scratch and reduced once, not once per edge. This amortizes the dominant edge
bandwidth across the whole prime batch -- the big win when the table is RAM\--resident.

@<Subroutines@>=
void spmv_run_batch(int Nto,u64* pr,int B){
  int cc,b; long i;
  if(c0_g<0 || Plevs<=0){
    for(cc=1;cc<=Nto;cc++) if(cc<=direct_valid_col_g){ printf("%d",cc);
      for(b=0;b<B;b++) printf(" %llu",(unsigned long long)(u64)(cnt[cc*m]%pr[b])); printf("\n"); }
    return; }
  u32* cnt2b=xcalloc((size_t)(1<<16)*B,sizeof(u32));
  u32* v=xmalloc((size_t)nstate[0]*B*sizeof(u32));
  for(i=0;i<nstate[0];i++){ u128 sd=seedv[i]; for(b=0;b<B;b++) v[i*B+b]=(u32)(sd%pr[b]); }
  int basecol=c0_g+1;
  static Edge* lbuf=0; static long lbufcap=0;
  while(basecol<=Nto){ int L;
    for(L=0;L<Plevs;L++){ int abscol=basecol+L/m, substep=L%m; long e;
      for(e=0;e<ncomp[L];e++){ int idx=abscol*m+substep+1+comps[L][e].delta; long sc=comps[L][e].src; u64 mult=comps[L][e].mult;
        if(idx>=0 && idx<(1<<16)) for(b=0;b<B;b++){ size_t o=(size_t)idx*B+b;
          cnt2b[o]=(u32)(((u64)cnt2b[o]+(u64)v[(size_t)sc*B+b]*mult)%pr[b]); } }
      u32* vn=xcalloc((size_t)nstate[L+1]*B,sizeof(u32));
      long ne=nedge[L];
      /* read the level's edges (served from page cache if resident), then apply
         in parallel: split at |dst| boundaries so threads write disjoint |vn| */
      if(lbufcap<ne){ lbufcap=ne; lbuf=xrealloc(lbuf,lbufcap*sizeof(Edge)); }
      fseeko(tbl_fp,lev_edge_off[L],SEEK_SET);
      { long rd=0; while(rd<ne){ long ch=ne-rd; if(ch>(1L<<24)) ch=1L<<24;
          if((long)fread(lbuf+rd,sizeof(Edge),ch,tbl_fp)!=ch){ fprintf(stderr,"table read error\n"); exit(3); } rd+=ch; } }
      { int NT=omp_get_max_threads(),tt; static long pbnd[NTMAX+1]; pbnd[0]=0; pbnd[NT]=ne;
        for(tt=1;tt<NT;tt++){ long g=(long)((double)tt*ne/NT);
          while(g>0&&g<ne&&lbuf[g].dst==lbuf[g-1].dst) g++; pbnd[tt]=g; }
#pragma omp parallel for schedule(dynamic,1)
        for(tt=0;tt<NT;tt++){ long e2,cur=-1; u64 acc[64]; int bb;
          for(e2=pbnd[tt];e2<pbnd[tt+1];e2++){ long dst=lbuf[e2].dst, src=lbuf[e2].src; u64 c=lbuf[e2].c;
            if(dst!=cur){ if(cur>=0) for(bb=0;bb<B;bb++) vn[(size_t)cur*B+bb]=(u32)(acc[bb]%pr[bb]);
              cur=dst; for(bb=0;bb<B;bb++) acc[bb]=0; }
            for(bb=0;bb<B;bb++) acc[bb]+=(u64)v[(size_t)src*B+bb]*c; }
          if(cur>=0) for(bb=0;bb<B;bb++) vn[(size_t)cur*B+bb]=(u32)(acc[bb]%pr[bb]); } }
      free(v); v=vn; }
    basecol+=period_g; }
  free(v);
  for(cc=1;cc<=Nto;cc++){ printf("%d",cc);
    for(b=0;b<B;b++){ u64 r = cc<=direct_valid_col_g ? (u64)(cnt[cc*m]%pr[b]) : (u64)cnt2b[(size_t)(cc*m)*B+b];
      printf(" %llu",(unsigned long long)r); } printf("\n"); }
  free(cnt2b); }

@ Compressed edge tables (scheme A: CSC-style delta + varint). Each level's
coalesced, |dst|-sorted edges are encoded as a byte stream of groups
$\langle$dst\--delta, gsize, (src\--delta, c)$^{gsize}\rangle$ (all LEB128 varints);
per\--group boundaries let us store a handful of split points (byte offset +
running dst) so the stream decodes in parallel. Typically $\sim3$ bytes/edge vs 16,
so a 100\,GB table shrinks to $\sim$25\,GB -- it then fits in RAM, turning the
disk\--bound SpMV into a cache\--bound one, and cuts memory traffic $\sim5\times$.

@<Globals@>=
unsigned char* cedge[MAXLEV]; long ced_len[MAXLEV];
long* csoff[MAXLEV]; long* csdst[MAXLEV]; int csn[MAXLEV];

@ @<Subroutines@>=
static int vput(unsigned char*p,u64 x){ int n=0; while(x>=0x80){ p[n++]=(x&0x7f)|0x80; x>>=7; } p[n++]=(unsigned char)x; return n; }
static int vlen(u64 x){ int n=1; while(x>=0x80){ x>>=7; n++; } return n; }   /* bytes vput would write */
static u64 vget(unsigned char**pp){ unsigned char*p=*pp; u64 x=0; int sh=0; unsigned char b; do{ b=*p++; x|=(u64)(b&0x7f)<<sh; sh+=7; }while(b&0x80); *pp=p; return x; }
static inline u32 barr(u64 x,u64 p,u64 mu){ u64 q=(u64)(((u128)x*mu)>>64); u64 r=x-q*p; if(r>=p)r-=p; if(r>=p)r-=p; return (u32)r; }
void compress_table(const char*in,const char*out){
  FILE*f=fopen(in,"rb"); if(!f){ fprintf(stderr,"cannot open %s\n",in); exit(1); }
  int hdr[6]; if(fread(hdr,sizeof(int),6,f)!=6){fprintf(stderr,"bad table\n");exit(1);}
  m=hdr[0];period_g=hdr[1];c0_g=hdr[2];Plevs=hdr[3];recend_col_g=hdr[4];direct_valid_col_g=hdr[5];
  fread(&seedn,sizeof(long),1,f); seedv=xmalloc(seedn*sizeof(u128)); fread(seedv,sizeof(u128),seedn,f);
  fread(nstate,sizeof(long),Plevs+1,f);
  FILE*g=fopen(out,"wb");
  fwrite(hdr,sizeof(int),6,g); fwrite(&seedn,sizeof(long),1,g); fwrite(seedv,sizeof(u128),seedn,g); fwrite(nstate,sizeof(long),Plevs+1,g);
  long SPLIT=1L<<21; Edge* eb=0; long ebcap=0; unsigned char* ob=0; long obcap=0; int L;
  for(L=0;L<Plevs;L++){
    long ne; fread(&ne,sizeof(long),1,f);
    if(ebcap<ne){ ebcap=ne; eb=xrealloc(eb,ebcap*sizeof(Edge)); }
    if(ne) fread(eb,sizeof(Edge),ne,f);
    if(obcap<ne*6+1024){ obcap=ne*6+1024; ob=xrealloc(ob,obcap); }
    long olen=0, prev_dst=-1, since=0, e2=0; int scap=(int)(ne/SPLIT)+2;
    long* soff=xmalloc(scap*sizeof(long)); long* sdst=xmalloc(scap*sizeof(long)); int sn=0;
    while(e2<ne){
      if(e2==0 || since>=SPLIT){ soff[sn]=olen; sdst[sn]=prev_dst; sn++; since=0; }
      long dst=eb[e2].dst; olen+=vput(ob+olen,(u64)(dst-prev_dst));
      long j=e2; while(j<ne && eb[j].dst==dst) j++; long gsize=j-e2;
      olen+=vput(ob+olen,(u64)gsize);
      long prev_src=0,k; for(k=e2;k<j;k++){ olen+=vput(ob+olen,(u64)(eb[k].src-prev_src)); prev_src=eb[k].src; olen+=vput(ob+olen,eb[k].c); }
      prev_dst=dst; since+=gsize; e2=j;
    }
    fwrite(&olen,sizeof(long),1,g); fwrite(&sn,sizeof(int),1,g);
    fwrite(soff,sizeof(long),sn,g); fwrite(sdst,sizeof(long),sn,g); if(olen) fwrite(ob,1,olen,g);
    free(soff); free(sdst);
    long ncp; fread(&ncp,sizeof(long),1,f); Comp*cp=xmalloc((ncp+1)*sizeof(Comp)); if(ncp) fread(cp,sizeof(Comp),ncp,f);
    fwrite(&ncp,sizeof(long),1,g); if(ncp) fwrite(cp,sizeof(Comp),ncp,g); free(cp);
  }
  int nc; fread(&nc,sizeof(int),1,f); memset(cnt,0,sizeof(cnt)); fread(cnt,sizeof(u128),nc,f);
  fwrite(&nc,sizeof(int),1,g); fwrite(cnt,sizeof(u128),nc,g);
  free(eb); free(ob); fclose(f); fclose(g);
  fprintf(stderr,"compressed %s -> %s\n",in,out);
}
int load_ctables(const char*path){ FILE*f=fopen(path,"rb"); if(!f) return 0; int L,hdr[6];
  if(fread(hdr,sizeof(int),6,f)!=6) return 0; m=hdr[0];period_g=hdr[1];c0_g=hdr[2];Plevs=hdr[3];recend_col_g=hdr[4];direct_valid_col_g=hdr[5];
  fread(&seedn,sizeof(long),1,f); seedv=xmalloc(seedn*sizeof(u128)); fread(seedv,sizeof(u128),seedn,f);
  fread(nstate,sizeof(long),Plevs+1,f);
  for(L=0;L<Plevs;L++){ fread(&ced_len[L],sizeof(long),1,f); fread(&csn[L],sizeof(int),1,f);
    csoff[L]=xmalloc((csn[L]+1)*sizeof(long)); csdst[L]=xmalloc((csn[L]+1)*sizeof(long));
    fread(csoff[L],sizeof(long),csn[L],f); fread(csdst[L],sizeof(long),csn[L],f);
    cedge[L]=xmalloc(ced_len[L]+16); fread(cedge[L],1,ced_len[L],f);
    fread(&ncomp[L],sizeof(long),1,f); comps[L]=xmalloc((ncomp[L]+1)*sizeof(Comp)); fread(comps[L],sizeof(Comp),ncomp[L],f); }
  int nc; fread(&nc,sizeof(int),1,f); memset(cnt,0,sizeof(cnt)); fread(cnt,sizeof(u128),nc,f);
  fclose(f); return 1; }
void spmv_run_batch_c(int Nto,u64* pr,int B){
  int cc,b; long i;
  if(c0_g<0||Plevs<=0){ for(cc=1;cc<=Nto;cc++) if(cc<=direct_valid_col_g){ printf("%d",cc);
    for(b=0;b<B;b++) printf(" %llu",(unsigned long long)(u64)(cnt[cc*m]%pr[b])); printf("\n"); } return; }
  u32* cnt2b=xcalloc((size_t)(1<<16)*B,sizeof(u32));
  u32* v=xmalloc((size_t)nstate[0]*B*sizeof(u32));
  for(i=0;i<nstate[0];i++){ u128 sd=seedv[i]; for(b=0;b<B;b++) v[i*B+b]=(u32)(sd%pr[b]); }
  static u64 mu[64]; for(b=0;b<B;b++) mu[b]=(u64)((((u128)1)<<64)/pr[b]);
  int basecol=c0_g+1;
  while(basecol<=Nto){ int L;
    for(L=0;L<Plevs;L++){ int abscol=basecol+L/m, substep=L%m; long e;
      for(e=0;e<ncomp[L];e++){ int idx=abscol*m+substep+1+comps[L][e].delta; long sc=comps[L][e].src; u64 mult=comps[L][e].mult;
        if(idx>=0&&idx<(1<<16)) for(b=0;b<B;b++){ size_t o=(size_t)idx*B+b; cnt2b[o]=(u32)(((u64)cnt2b[o]+(u64)v[(size_t)sc*B+b]*mult)%pr[b]); } }
      u32* vn=xcalloc((size_t)nstate[L+1]*B,sizeof(u32));
      int NT=omp_get_max_threads(),tt, SN=csn[L];
#pragma omp parallel for schedule(dynamic,1)
      for(tt=0;tt<NT;tt++){ int s0=(int)((long)tt*SN/NT), s1=(int)((long)(tt+1)*SN/NT); if(s1>SN)s1=SN;
        if(s0<s1){ unsigned char* p=cedge[L]+csoff[L][s0];
          unsigned char* endp=cedge[L]+(s1<SN? csoff[L][s1] : ced_len[L]);
          long prev_dst=csdst[L][s0], cur=-1; u64 acc[64]; int bb,w;
          long Wdst[640],Wsrc[640]; u32 Wc[640];
          while(p<endp){ int nw=0;               /* decode a window of whole groups */
            while(p<endp && nw<512){ long dst=prev_dst+(long)vget(&p); long gsize=(long)vget(&p); long ps=0,gg;
              for(gg=0;gg<gsize;gg++){ long src=ps+(long)vget(&p); ps=src; Wdst[nw]=dst; Wsrc[nw]=src; Wc[nw]=(u32)vget(&p); nw++; }
              prev_dst=dst; }
            for(w=0;w<nw;w++) __builtin_prefetch(&v[(size_t)Wsrc[w]*B],0,1);  /* hide gather latency */
            for(w=0;w<nw;w++){ if(Wdst[w]!=cur){ if(cur>=0) for(bb=0;bb<B;bb++) vn[(size_t)cur*B+bb]=barr(acc[bb],pr[bb],mu[bb]);
                cur=Wdst[w]; for(bb=0;bb<B;bb++) acc[bb]=0; }
              for(bb=0;bb<B;bb++) acc[bb]+=(u64)v[(size_t)Wsrc[w]*B+bb]*Wc[w]; } }
          if(cur>=0) for(bb=0;bb<B;bb++) vn[(size_t)cur*B+bb]=(u32)(acc[bb]%pr[bb]); } }
      free(v); v=vn; }
    basecol+=period_g; }
  free(v);
  for(cc=1;cc<=Nto;cc++){ printf("%d",cc);
    for(b=0;b<B;b++){ u64 r=cc<=direct_valid_col_g?(u64)(cnt[cc*m]%pr[b]):(u64)cnt2b[(size_t)(cc*m)*B+b]; printf(" %llu",(unsigned long long)r); } printf("\n"); }
  free(cnt2b); }

@* Main.
No arguments: self\--checks (direct sweep vs known values, and SpMV vs direct).
|dualham m Wb N|: build to column $W_b$, then extend to column~$N$ by SpMV.

@<Main@>=
void spmv_run(int Nto,int crosscheck);
void spmv_run_batch(int Nto,u64* pr,int B);
void compress_table(const char*,const char*);
int load_ctables(const char*);
void spmv_run_batch_c(int Nto,u64* pr,int B);
int build_extract(int mm,int Wb);
@* Exact big integers in one shot: build once, CRT over hard\--wired primes.
The build extracts the (prime\--agnostic) edge tables once; we then replay the
weight\--free SpMV once per prime |p_b| (mod |p_b|) and Garner\--combine the
residues into an exact big integer per column -- all in RAM, no disk. The number
of primes is fixed by a rigorous bound: an open tour is a permutation of the
$mn$ cells whose consecutive cells are knight\--adjacent (degree $\le 8$), so
there are fewer than $8^{mn}<2^{3mn}$ of them; with 31\--bit primes that needs
$K=\lceil 3mn/31\rceil$, and we take a few extra so that {\it dropping\/} the
extras leaves every value unchanged (an empirical confirmation of exactness).

@<Subroutines@>=
static u64 modpowu(u64 a,u64 e,u64 p){ u64 x=1; a%=p; while(e){ if(e&1) x=(u128)x*a%p; a=(u128)a*a%p; e>>=1; } return x; }
static int isprime(u64 n){ u64 q; if(n<2) return 0; for(q=2;q<40&&q*q<=n;q++) if(n%q==0) return n==q;
  u64 d=n-1; int r=0,i,j; while(!(d&1)){ d>>=1; r++; }
  u64 A[]={2,3,5,7,11,13,17,19,23,29,31,37};
  for(i=0;i<12;i++){ u64 a=A[i]%n; if(!a) continue; u64 x=modpowu(a,d,n); if(x==1||x==n-1) continue; int ok=0;
    for(j=0;j<r-1;j++){ x=(u128)x*x%n; if(x==n-1){ok=1;break;} } if(!ok) return 0; }
  return 1; }
static int gen_primes(u64*pr,int K){ int n=0; u64 x=(1ULL<<31)-1; while(n<K && x>2){ if(isprime(x)) pr[n++]=x; x-=2; } return n; }
static u64 modinv(u64 a,u64 p){ return modpowu(a%p,p-2,p); }
@#
#define BNCAP 4096
typedef struct{ u32 d[BNCAP]; int n; } BN;   /* little-endian, base $10^9$ */
static void bn_setu64(BN*x,u64 v){ x->n=0; while(v){ x->d[x->n++]=(u32)(v%1000000000ULL); v/=1000000000ULL; } if(!x->n){x->d[0]=0;x->n=1;} }
static u32 bn_mod(const BN*x,u64 p){ u64 r=0; int i; for(i=x->n-1;i>=0;i--) r=(r*1000000000ULL + x->d[i])%p; return (u32)r; }
static void bn_mulsmall(BN*x,u64 s){ u64 c=0; int i; for(i=0;i<x->n;i++){ u64 v=(u64)x->d[i]*s+c; x->d[i]=(u32)(v%1000000000ULL); c=v/1000000000ULL; } while(c){ x->d[x->n++]=(u32)(c%1000000000ULL); c/=1000000000ULL; } }
static void bn_addmul(BN*x,const BN*M,u64 t){ u64 c=0; int i,n=x->n>M->n?x->n:M->n; for(i=0;i<n||c;i++){ u64 v=c; if(i<x->n)v+=x->d[i]; if(i<M->n)v+=(u64)M->d[i]*t; x->d[i]=(u32)(v%1000000000ULL); c=v/1000000000ULL; } if(i>x->n)x->n=i; while(x->n>1 && x->d[x->n-1]==0) x->n--; }
static void bn_print(const BN*x){ int i; printf("%u",x->d[x->n-1]); for(i=x->n-2;i>=0;i--) printf("%09u",x->d[i]); }

@ Batched weight-free replay: apply |B| primes in ONE pass over the (compressed)
edge tables, with interleaved |u32| vectors |v[i*B+b]| and Barrett reduction --
so the expensive edge decode is amortised across |B| residues (|K/B| passes for
|K| primes instead of |K|). Fills |cnt2b_g| (residues, indexed |(col)*B+b|).

@ @<Subroutines@>=
/* Fills |out| (caller-owned, size |vcols*B|) with residues indexed |col*B+b|;
   uses no globals, so a batch of primes runs on its own thread. */
void wfree_replay_batch(int Nto,u64*pr,int B,u32*out,size_t vcols){
  int b; u64 mu[64]; for(b=0;b<B;b++) mu[b]=(u64)((((u128)1)<<64)/pr[b]);
  memset(out,0,vcols*B*sizeof(u32));
  u32* cnt2b_g=out;
  u32* v=xcalloc((size_t)nstate[0]*B,sizeof(u32)); for(b=0;b<B;b++) v[b]=1;
/* Parallel decode of one level: the split points |csoff|/|csdst| cut the byte
   stream at dst-group boundaries into disjoint dst blocks, so each thread decodes
   its span and scatters into |vn_| with no collisions (window + prefetch + Barrett). */
#define APPLYB(LL,ACOL,SUB) do{ int L_=(LL),abscol_=(ACOL),substep_=(SUB); long e_; \
    for(e_=0;e_<ncomp[L_];e_++){ size_t idx_=(size_t)(abscol_*m+substep_+1+comps[L_][e_].delta)*B; long sc_=comps[L_][e_].src; u64 mult_=comps[L_][e_].mult; \
      for(b=0;b<B;b++){ cnt2b_g[idx_+b]=barr((u64)cnt2b_g[idx_+b]+(u64)v[(size_t)sc_*B+b]*mult_,pr[b],mu[b]); } } \
    u32* vn_=xcalloc((size_t)nstate[L_+1]*B,sizeof(u32)); \
    int NT_=omp_get_max_threads(), tt_, SN_=csn[L_]; \
    _Pragma("omp parallel for schedule(dynamic,1)") \
    for(tt_=0;tt_<NT_;tt_++){ int s0_=(int)((long)tt_*SN_/NT_), s1_=(int)((long)(tt_+1)*SN_/NT_); if(s1_>SN_)s1_=SN_; \
      if(s0_<s1_){ unsigned char* p_=cedge[L_]+csoff[L_][s0_], *end_=cedge[L_]+(s1_<SN_?csoff[L_][s1_]:ced_len[L_]); \
        long prev_=csdst[L_][s0_], cur_=-1; u64 acc_[64]; int bb_,w_; long Wd_[640],Ws_[640]; u32 Wc_[640]; \
        while(p_<end_){ int nw_=0; \
          while(p_<end_ && nw_<512){ long dst=prev_+(long)vget(&p_); long gs=(long)vget(&p_),ps=0,gg; \
            for(gg=0;gg<gs;gg++){ long src=ps+(long)vget(&p_); ps=src; Wd_[nw_]=dst; Ws_[nw_]=src; Wc_[nw_]=(u32)vget(&p_); nw_++; } prev_=dst; } \
          for(w_=0;w_<nw_;w_++) __builtin_prefetch(&v[(size_t)Ws_[w_]*B],0,1); \
          for(w_=0;w_<nw_;w_++){ if(Wd_[w_]!=cur_){ if(cur_>=0) for(bb_=0;bb_<B;bb_++) vn_[(size_t)cur_*B+bb_]=barr(acc_[bb_],pr[bb_],mu[bb_]); cur_=Wd_[w_]; for(bb_=0;bb_<B;bb_++) acc_[bb_]=0; } \
            for(bb_=0;bb_<B;bb_++) acc_[bb_]+=(u64)v[(size_t)Ws_[w_]*B+bb_]*Wc_[w_]; } } \
        if(cur_>=0) for(bb_=0;bb_<B;bb_++) vn_[(size_t)cur_*B+bb_]=(u32)(acc_[bb_]%pr[bb_]); } } \
    free(v); v=vn_; }while(0)
  int L, X, ss, pfx=(aggr_R>=0)?junc_P0:Ltot_g, perbase=(aggr_R>=0)?junc_P0:(Ltot_g-Plevs), pref=c0_g+1;
  for(L=0;L<pfx;L++) APPLYB(L,L/m,L%m);
  int basecol;
  if(aggr_R>=0){  /* remap the col-R boundary vector (col-R order) into the stable order, then run the periodic block from col R+1 */
    u32* vs=xcalloc((size_t)junc_nstab*B,sizeof(u32)); long i2;
    for(i2=0;i2<junc_nkr;i2++){ int j=junc_remap[i2]; if(j>=0){ int bb; for(bb=0;bb<B;bb++) vs[(size_t)j*B+bb]=v[(size_t)i2*B+bb]; } }
    free(v); v=vs; basecol=aggr_R+1;
  } else basecol=recend_col_g+1;
  for(X=basecol; X<=Nto; X++){ int pc=((X-pref)%period_g+period_g)%period_g;
    for(ss=0;ss<m;ss++) APPLYB(perbase+pc*m+ss, X, ss); }
  free(v);
#undef APPLYB
}

@ @<Subroutines@>=
static double nowsec(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+t.tv_nsec*1e-9; }
void exact_mode(int mm,int N){
  int Wb=20, b, c; int K=(3*mm*N)/31 + 6;   /* rigorous bound + slack for the drop-a-prime check */
  static u64 pr[512]; if(K>500)K=500; int nP=gen_primes(pr,K);
  fprintf(stderr,"exact m=%d to n=%d : %d primes (product > 2^%d, bound < 2^%d)\n",mm,N,nP,31*nP,3*mm*N);
  double t0=nowsec();
  no_replay=1; MODP=0; run_periodic(mm,Wb,N);   /* build the edge tables in RAM, no replay */
  double tb=nowsec(); fprintf(stderr,"[time] build %.1fs  (expand %.1fs, reduce %.1fs [of which compress %.1fs])\n",
    tb-t0, t_expand, t_reduce, t_compress);
  no_replay=0;
  static u32* resid=0; resid=xrealloc(resid,(size_t)nP*(N+1)*sizeof(u32));
  long maxns=1; { int L; for(L=0;L<=Ltot_g;L++) if(nstate[L]>maxns) maxns=nstate[L]; }
  size_t vcols=(size_t)((N+3)*m+8);
  /* The per-level decode is now PARALLEL (split points), so each batch already uses
     all cores. Make Bp as LARGE as memory allows (max decode amortisation) and run
     the few batches SERIALLY -- the decode fills the cores within each batch. */
  long cap = mem_cap_bytes>0? mem_cap_bytes : (100L<<30);
  long budget = cap - rss_bytes() - (8L<<30);   /* leave 8GB headroom */
  int Bp=(int)(budget/((long)maxns*8L+(1L<<20)));   /* v+vn interleaved: 2 x maxns x Bp x 4B */
  if(Bp<1)Bp=1; if(Bp>nP)Bp=nP; if(Bp>64)Bp=64;   /* Barrett acc[64] limits B<=64 */
  int nbatch=(nP+Bp-1)/Bp, bi;
  fprintf(stderr,"replay: %d primes, %d batches of %d, parallel decode (maxns=%ld, budget %.0fGB)\n",nP,nbatch,Bp,maxns,budget/1073741824.0);
  for(bi=0;bi<nbatch;bi++){ int base=bi*Bp, bb=nP-base; if(bb>Bp)bb=Bp;
    u32* out=xcalloc(vcols*bb,sizeof(u32));
    wfree_replay_batch(N,pr+base,bb,out,vcols);
    int cc,k; for(cc=1;cc<=N;cc++) for(k=0;k<bb;k++) resid[(size_t)(base+k)*(N+1)+cc]=out[(size_t)(cc*m)*bb+k];
    free(out); }
  { double tr=nowsec(); fprintf(stderr,"[time] replay %.1fs\n",tr-tb); tb=tr; }
  for(c=1;c<=N;c++){
    BN x,M; bn_setu64(&x,resid[c]); bn_setu64(&M,pr[0]);
    for(b=1;b<nP;b++){ u64 p=pr[b], xm=bn_mod(&x,p);
      u64 t=(u64)((u128)((resid[(size_t)b*(N+1)+c]+p-xm)%p) * modinv(bn_mod(&M,p),p) % p);
      bn_addmul(&x,&M,t); bn_mulsmall(&M,p); }
    printf("open %dx%d = ",mm,c); bn_print(&x); printf("\n");
  }
}

@ @<Subroutines@>=
int main(int argc,char*argv[]){
  start_mem_watcher();
  if(argc==4 && !strcmp(argv[1],"exact")){ exact_mode(atoi(argv[2]),atoi(argv[3])); return 0; }  /* exact m N : one program, from scratch to open m x 1..N as exact big integers, all in RAM */
  if(argc>=5 && !strcmp(argv[1],"build")){ /* build m Wb file [ckpt] : build+extract, dump tables */
    stop_after_record=1; stream_edges=1; snprintf(rec_path,sizeof rec_path,"%s",argv[4]);
    if(argc>=6) ckpt_path=argv[5]; if(argc>=7) ckpt_every=atoi(argv[6]);
    build_extract(atoi(argv[2]),atoi(argv[3])); dump_tables(argv[4]); return 0; }
  if(argc==5 && !strcmp(argv[1],"resume")){ /* resume ckpt Wb file : continue a build */
    stop_after_record=1; stream_edges=1; snprintf(rec_path,sizeof rec_path,"%s",argv[4]);
    resume_from=load_ckpt(argv[2]); ckpt_path=argv[2];
    build_extract(m,atoi(argv[3])); dump_tables(argv[4]); return 0; }
  if(argc>=4 && !strcmp(argv[1],"run")){ /* run file Nto [prime] : load tables, SpMV (mod prime if given) */
    if(!load_tables(argv[2])){ fprintf(stderr,"cannot load %s\n",argv[2]); return 1; }
    if(argc>=5) MODP=strtoull(argv[4],0,10);
    spmv_run(atoi(argv[3]),0); return 0; }
  if(argc>=5 && !strcmp(argv[1],"runb")){ /* runb file Nto p1 p2 ... : one table pass, many primes */
    if(!load_tables(argv[2])){ fprintf(stderr,"cannot load %s\n",argv[2]); return 1; }
    int B=argc-4; if(B>64) B=64; static u64 pr[64]; int b; for(b=0;b<B;b++) pr[b]=strtoull(argv[4+b],0,10);
    spmv_run_batch(atoi(argv[3]),pr,B); return 0; }
  if(argc==4 && !strcmp(argv[1],"compress")){ compress_table(argv[2],argv[3]); return 0; }
  if(argc>=5 && !strcmp(argv[1],"runbc")){ /* runbc cfile Nto p1 p2 ... : compressed table, batched */
    if(!load_ctables(argv[2])){ fprintf(stderr,"cannot load %s\n",argv[2]); return 1; }
    int B=argc-4; if(B>64) B=64; static u64 pr[64]; int b; for(b=0;b<B;b++) pr[b]=strtoull(argv[4+b],0,10);
    spmv_run_batch_c(atoi(argv[3]),pr,B); return 0; }
  if(argc==4){ run_periodic(atoi(argv[1]),atoi(argv[2]),atoi(argv[3])); return 0; }
  /* A light smoke test: m=4 is tiny, and by transposition it also pins the m=6
     and m=7 values (4xc=cx4) without their million-state sweeps; one m=5 strip,
     built wide enough to record, exercises both the direct and the SpMV paths. */
  struct{int m,Wb,Next,c;u64 e;} chk[]={
    {4,12,9,6,744ULL},                      /* 4x6 = 6x4 */
    {4,12,9,7,6378ULL},                     /* 4x7 = 7x4 */
    {5,18,22,8,18061054ULL},                /* m=5, direct column */
    {5,18,22,12,3611823644006ULL},          /* m=5, direct seam column (recend+1) */
    {5,18,22,15,24535910156176100ULL},      /* m=5, SpMV extrapolation (past recend+1) */
    {0,0,0,0,0}};
  int i,bad=0,lm=0,lw=0,ln=0;
  for(i=0;chk[i].m;i++){
    if(chk[i].m!=lm||chk[i].Wb!=lw||chk[i].Next!=ln){
      run_periodic(chk[i].m,chk[i].Wb,chk[i].Next); lm=chk[i].m; lw=chk[i].Wb; ln=chk[i].Next; }
    int c=chk[i].c; u64 g=cnt2[c*m], e=red(chk[i].e);   /* weight-free: all counts come from the replay */
    printf("open %dx%d = %llu  exp %llu  %s\n",chk[i].m,c,g,e,g==e?"OK":"FAIL");
    if(g!=e) bad++; }
  printf("%s\n",bad?"SOME FAILED":"ALL OK"); return bad?1:0;
}
