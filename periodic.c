/* periodic.c -- extract the PERIODIC transfer of the dual-frontier engine, then
   compute open m x c by a sparse matrix-vector product (SpMV) instead of
   re-running generate every column.  Validates against the direct sweep.
   Because the dual key is translation-invariant, the sorted bucket at each
   period-substep is identical every period, so bucket indices are stable ids
   and the extracted integer edge tables are reusable. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define MAXF 512
#define MAXLEV 40

int m,n,V;
int NB[64*64][8], ND[64*64];
static const int KR[8]={-2,-2,-1,-1,1,1,2,2};
static const int KC[8]={-1,1,-2,2,-2,2,-1,1};
void build(void){ int r,c,k; V=m*n;
  for(c=0;c<n;c++)for(r=0;r<m;r++){ int v=c*m+r; ND[v]=0;
    for(k=0;k<8;k++){ int rr=r+KR[k],cc=c+KC[k];
      if(rr>=0&&rr<m&&cc>=0&&cc<n) NB[v][ND[v]++]=cc*m+rr; } } }
int fr[MAXF], ifrb[64*64+2];
int frontier_before(int s){ int q=0,v,u,k; static char inF[64*64+2];
  for(v=0;v<=V;v++) inF[v]=0; inF[V]=1; if(s<V) inF[s]=1;
  for(v=s+1;v<V;v++) for(k=0;k<ND[v];k++){ u=NB[v][k]; if(u<s){ inF[v]=1; break; } }
  for(v=0;v<=V;v++) if(inF[v]){ fr[q]=v; ifrb[v]=q; q++; }
  return q; }
int mate[MAXF]; int cycle;
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

/* bucket */
unsigned char*kp; long kcap,kuse; typedef struct{long off;int len;u64 w; long src;} Rec; Rec*rc; long nr,rcap;
unsigned char*cb;
int rcmp(const void*A,const void*B){const Rec*a=A,*b=B;int l=a->len<b->len?a->len:b->len;int d=memcmp(cb+a->off,cb+b->off,l);return d?d:a->len-b->len;}
void emit(unsigned char*key,int len,u64 w,long src){
  if(kuse+len>kcap){kcap=kcap*2+len+65536;kp=realloc(kp,kcap);}
  if(nr>=rcap){rcap=rcap*2+65536;rc=realloc(rc,rcap*sizeof(Rec));}
  memcpy(kp+kuse,key,len);rc[nr].off=kuse;rc[nr].len=len;rc[nr].w=w;rc[nr].src=src;kuse+=len;nr++; }
unsigned char*curkp; long*curoff; int*curkl; u64*curw; long ncur;
u64 cnt[1<<16];
int qnew,posS,apexnew,STEMP; static int bmate[MAXF],o2n[MAXF];
void build_bmate(int*omate,int qold){ int i;
  for(i=0;i<=qnew;i++) bmate[i]=-2; STEMP=qnew;
  for(i=0;i<qold;i++){ int dst=(i==posS)?STEMP:o2n[i]; if(dst<0) continue;
    if(omate[i]==-1) bmate[dst]=-1;
    else if(omate[i]>=0){ int op=omate[i]; int pdst=(op==posS)?STEMP:o2n[op]; bmate[dst]=pdst>=0?pdst:-2; } } }

/* ---- periodic recording ---- */
int recording=0, reclev=0;
typedef struct{ int src,dst; u64 c; } Edge;
Edge* edges[MAXLEV]; long nedge[MAXLEV];
typedef struct{ int src,delta; u64 mult; } Comp;
Comp* comps[MAXLEV]; long ncomp[MAXLEV];
long nstate[MAXLEV+1];        /* #states at each recorded level */
u64* seedv; long seedn;
int Plevs;                   /* total recorded levels = period*m */
/* scratch during a recorded step */
Comp reccomp_buf[1<<20]; long reccomp_n;

int ecmp(const void*A,const void*B){const Edge*a=A,*b=B; if(a->src!=b->src)return a->src-b->src; return a->dst-b->dst;}

void run_step(int s,int rec){
  int i; long si; long in_ncur=ncur; int qold=frontier_before(s); posS=ifrb[s];
  static int frold[MAXF]; for(i=0;i<qold;i++) frold[i]=fr[i];
  qnew=frontier_before(s+1); static int ifrnew[64*64+2], frnew[MAXF];
  for(i=0;i<qnew;i++) frnew[i]=fr[i]; for(i=0;i<=V;i++) ifrnew[i]=-1; for(i=0;i<qnew;i++) ifrnew[frnew[i]]=i;
  apexnew=ifrnew[V];
  int nbr[16],rr=0; nbr[rr++]=apexnew;
  for(i=0;i<ND[s];i++){ int w=NB[s][i]; if(w>s && ifrnew[w]>=0) nbr[rr++]=ifrnew[w]; }
  for(i=0;i<qold;i++) o2n[i]=(frold[i]==s)?-1:ifrnew[frold[i]];
  nr=0; kuse=0; reccomp_n=0;
  for(si=0;si<ncur;si++){
    unsigned char*ok=curkp+curoff[si]; int okl=curkl[si]; u64 w=curw[si];
    static int omate[MAXF]; int a,b;
    for(i=0;i<okl;i++){int c=ok[i]; omate[i]= c==0?-2 : c==255?-1 : c-1;}
    int deg=omate[posS]==-2?0:omate[posS]==-1?2:1, need=2-deg;
    unsigned char nk[MAXF]; int nl;
    if(need==0){ build_bmate(omate,qold); for(i=0;i<qnew;i++) mate[i]=bmate[i]; nl=keyof(qnew,nk); emit(nk,nl,w,si); }
    else if(need==1){ for(a=0;a<rr;a++){ build_bmate(omate,qold); for(i=0;i<=qnew;i++) mate[i]=bmate[i];
      if(add_derived(STEMP,nbr[a])){ nl=keyof(qnew,nk); emit(nk,nl,w,si); }
      else if(cycle){ int mp=completion_mp(qnew,apexnew,frnew,s); if(mp){ cnt[mp]+=w;
        if(rec){ reccomp_buf[reccomp_n].src=si; reccomp_buf[reccomp_n].delta=mp-(s+1); reccomp_buf[reccomp_n].mult=1; reccomp_n++; } } } } }
    else { for(a=0;a<rr;a++)for(b=a+1;b<rr;b++){ build_bmate(omate,qold); for(i=0;i<=qnew;i++) mate[i]=bmate[i];
      if(!add_derived(STEMP,nbr[a])){ if(cycle){ int mp=completion_mp(qnew,apexnew,frnew,s); if(mp){ cnt[mp]+=w;
        if(rec){ reccomp_buf[reccomp_n].src=si; reccomp_buf[reccomp_n].delta=mp-(s+1); reccomp_buf[reccomp_n].mult=1; reccomp_n++; } } } continue; }
      if(add_derived(STEMP,nbr[b])){ nl=keyof(qnew,nk); emit(nk,nl,w,si); }
      else if(cycle){ int mp=completion_mp(qnew,apexnew,frnew,s); if(mp){ cnt[mp]+=w;
        if(rec){ reccomp_buf[reccomp_n].src=si; reccomp_buf[reccomp_n].delta=mp-(s+1); reccomp_buf[reccomp_n].mult=1; reccomp_n++; } } } } }
  }
  /* sort+reduce; if recording, capture edges (src bucket-index -> dst bucket-index) */
  cb=kp; qsort(rc,nr,sizeof(Rec),rcmp);
  curoff=realloc(curoff,(nr+1)*sizeof(long)); curkl=realloc(curkl,(nr+1)*sizeof(int));
  curw=realloc(curw,(nr+1)*sizeof(u64)); curkp=realloc(curkp,kuse+1);
  Edge* eb=0; long enb=0;
  if(rec){ eb=malloc(nr*sizeof(Edge)); }
  long k=0,use=0,ri;
  for(ri=0;ri<nr;){ long j=ri+1; u64 sw=rc[ri].w;
    if(rec){ long t; for(t=ri;t<nr;t++){ if(t>ri && !(rc[t].len==rc[ri].len && memcmp(kp+rc[t].off,kp+rc[ri].off,rc[ri].len)==0)) break;
        eb[enb].src=(int)rc[t].src; eb[enb].dst=(int)k; eb[enb].c=1; enb++; } }
    while(j<nr&&rc[j].len==rc[ri].len&&memcmp(kp+rc[j].off,kp+rc[ri].off,rc[ri].len)==0){sw+=rc[j].w;j++;}
    memcpy(curkp+use,kp+rc[ri].off,rc[ri].len); curoff[k]=use; curkl[k]=rc[ri].len; curw[k]=sw; use+=rc[ri].len; k++; ri=j; }
  ncur=k;
  if(rec){
    /* coalesce edges by (src,dst) */
    qsort(eb,enb,sizeof(Edge),ecmp);
    long o=0,p; for(p=0;p<enb;){ long q2=p+1; u64 cc=1;
      while(q2<enb&&eb[q2].src==eb[p].src&&eb[q2].dst==eb[p].dst){cc++;q2++;}
      eb[o]=eb[p]; eb[o].c=cc; o++; p=q2; }
    edges[reclev]=realloc(eb,o*sizeof(Edge)); nedge[reclev]=o;
    nstate[reclev]=in_ncur;      /* input (src) level size */
    nstate[reclev+1]=ncur;       /* output (dst) level size */
    /* comps */
    comps[reclev]=malloc((reccomp_n+1)*sizeof(Comp)); memcpy(comps[reclev],reccomp_buf,reccomp_n*sizeof(Comp)); ncomp[reclev]=reccomp_n;
    reclev++;
  }
}

int main(int argc,char*argv[]){
  if(argc<4){fprintf(stderr,"usage: periodic m Wbuild Nextend\n");return 1;}
  m=atoi(argv[1]); int Wb=atoi(argv[2]), Next=atoi(argv[3]); n=Wb; build();
  int i,s; int q=frontier_before(0); for(i=0;i<q;i++) mate[i]=-2;
  unsigned char key[MAXF]; int kl=keyof(q,key);
  curkp=malloc(1<<20); curoff=malloc(sizeof(long)); curkl=malloc(sizeof(int)); curw=malloc(sizeof(u64));
  memcpy(curkp,key,kl); curoff[0]=0; curkl[0]=kl; curw[0]=1; ncur=1;
  /* period detection by boundary key-set fingerprint */
  static u64 colfp[4096]; int period=0,c0=-1;
  int recstart=-1, recend=-1;
  for(s=0;s<V;s++){
    int atboundary = (s%m==m-1);
    int doing_rec = recording && s>=recstart && s<=recend;
    if(recording && s==recstart){ /* capture seed = current bucket weights (level 0) */
      seedn=ncur; seedv=malloc(ncur*sizeof(u64)); for(i=0;i<ncur;i++) seedv[i]=curw[i]; nstate[0]=ncur; reclev=0; }
    run_step(s, doing_rec);
    if(atboundary){
      int c=s/m; u64 h=1469598103934665603ULL; long t;
      for(t=0;t<ncur;t++){ unsigned char*kk=curkp+curoff[t]; int L=curkl[t],z; for(z=0;z<L;z++){ h^=kk[z]; h*=1099511628211ULL; } h^=0x9e; h*=1099511628211ULL; }
      colfp[c]=h;
      if(c0<0){ if(c>=1 && colfp[c]==colfp[c-1]){period=1;c0=c;} else if(c>=2 && colfp[c]==colfp[c-2]){period=2;c0=c;}
        if(c0>=0){ recstart=(c0+1)*m; recend=(c0+1+period)*m-1; Plevs=period*m;
          fprintf(stderr,"stable col %d period %d; recording cells %d..%d (%d levels)\n",c0,period,recstart,recend,Plevs);
          recording=1; }
      }
    }
  }
  /* reference open m x c from the direct sweep */
  /* SpMV from seed using recorded tables */
  static u64 cnt2[1<<16];
  u64* v=malloc(nstate[0]*sizeof(u64)); memcpy(v,seedv,nstate[0]*sizeof(u64));
  int basecol=c0+1; /* seed is boundary of col c0; first recorded level advances into col c0+1 */
  while(basecol<=Next){
    int L; for(L=0;L<Plevs;L++){
      int abscol=basecol + L/m, substep=L%m; long e;
      for(e=0;e<ncomp[L];e++){ int idx=abscol*m+substep+1+comps[L][e].delta; cnt2[idx]+=v[comps[L][e].src]*comps[L][e].mult; }
      u64* vn=calloc(nstate[L+1],sizeof(u64));
      for(e=0;e<nedge[L];e++) vn[edges[L][e].dst]+=v[edges[L][e].src]*edges[L][e].c;
      free(v); v=vn;
    }
    basecol+=period;
  }
  /* compare cnt (direct) vs cnt2 (SpMV) where the SpMV is fully covered.
     Completions crediting cnt[c*m] can occur up to ~m cells before c*m, so the
     SpMV (which starts recording at col c0+1) is complete from c0+3 onward. */
  int good=1,c, cfrom=c0+3;
  for(c=cfrom;c<=Next && c<=Wb;c++){
    if(cnt[c*m]!=cnt2[c*m]){ printf("MISMATCH c=%d direct=%llu spmv=%llu\n",c,cnt[c*m],cnt2[c*m]); good=0; }
    else printf("check c=%d: %llu  (direct==SpMV)\n",c,cnt[c*m]);
  }
  for(c=1;c<=Next;c++){ u64 val = c<=Wb? cnt[c*m] : cnt2[c*m]; if(val) printf("open %dx%d = %llu%s\n",m,c,val, c>Wb?"  (SpMV)":""); }
  printf("%s\n", good?"SpMV matches direct":"SpMV MISMATCH");
  return good?0:1;
}
