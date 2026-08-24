/* dual.c -- DUAL-frontier bucket transfer for open m x n knight tours.
   Frontier = UNPROCESSED cells adjacent to processed, plus an apex (open paths
   <-> cycles through the apex).  State = mate[] over the frontier (0 bare /
   -1 inner / else partner index).  DYNAHAM's compact representation, bucketed.
   Processing vertex s: give s a temporary slot, then add its (2-deg) config
   edges to future neighbours/apex via add_derived, so s ends at degree 2. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define MAXF 512

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
  for(v=0;v<=V;v++) inF[v]=0;
  inF[V]=1; if(s<V) inF[s]=1;
  for(v=s+1;v<V;v++) for(k=0;k<ND[v];k++){ u=NB[v][k]; if(u<s){ inF[v]=1; break; } }
  for(v=0;v<=V;v++) if(inF[v]){ fr[q]=v; ifrb[v]=q; q++; }
  return q; }

int mate[MAXF]; int cycle;
int add_derived(int i,int j){          /* bare=-2, inner=-1, outer=partner pos */
  cycle=0;
  if(mate[i]==-1||mate[j]==-1) return 0;
  if(mate[i]==-2){
    if(mate[j]==-2){ if(i==j){cycle=1;return 0;} mate[i]=j;mate[j]=i;return 1; }
    mate[i]=mate[j]; mate[mate[j]]=i; mate[j]=-1; return 1;
  } else if(mate[j]==-2){ mate[j]=mate[i]; mate[mate[i]]=j; mate[i]=-1; return 1; }
  else if(mate[i]!=j){ mate[mate[i]]=mate[j]; mate[mate[j]]=mate[i]; mate[i]=mate[j]=-1; return 1; }
  mate[i]=mate[j]=-1; cycle=1; return 0; }  /* joined ends of one subpath: cycle */

int keyof(int q,unsigned char*key){ int i;
  for(i=0;i<q;i++) key[i]= mate[i]==-2?0 : mate[i]==-1?255 : (unsigned char)(1+mate[i]);
  return q; }

/* cycle just closed: valid completion iff the apex is inner (in the loop) and no
   OUTER cells remain (every subpath consumed into the one cycle). Inner frontier
   cells are fine (covered, still awaiting their unprocessed neighbours). */
/* A cycle through the apex has closed.  It is a valid Hamiltonian m'-path iff
   the covered (inner) cells are exactly the contiguous prefix {0..m'-1} and every
   other board frontier cell is bare.  Returns m' (>0) to credit, else 0. */
int completion_mp(int q,int apexpos,int*frn,int s){ int k,mp=s+1;
  if(mate[apexpos]!=-1) return 0;
  k=0;
  while(k<q && frn[k]!=V && mate[k]==-1 && frn[k]==mp){ mp++; k++; }
  for(;k<q;k++){ if(frn[k]==V) continue; if(mate[k]!=-2) return 0; } /* rest bare */
  return mp;
}

unsigned char*kp; long kcap,kuse; typedef struct{long off;int len;u64 w;}Rec; Rec*rc; long nr,rcap;
unsigned char*cb;
int rcmp(const void*A,const void*B){const Rec*a=A,*b=B;int l=a->len<b->len?a->len:b->len;int d=memcmp(cb+a->off,cb+b->off,l);return d?d:a->len-b->len;}
void emit(unsigned char*key,int len,u64 w){
  if(kuse+len>kcap){kcap=kcap*2+len+65536;kp=realloc(kp,kcap);}
  if(nr>=rcap){rcap=rcap*2+65536;rc=realloc(rc,rcap*sizeof(Rec));}
  memcpy(kp+kuse,key,len);rc[nr].off=kuse;rc[nr].len=len;rc[nr].w=w;kuse+=len;nr++; }

unsigned char*curkp; long*curoff; int*curkl; u64*curw; long ncur;
u64 cnt[8192];

/* build base bmate over qnew real positions + one temp slot (=qnew) for s */
int qnew; static int bmate[MAXF]; static int o2n[MAXF]; int posS,apexnew,STEMP;
void build_bmate(int*omate,int qold){ int i;
  for(i=0;i<=qnew;i++) bmate[i]=-2;          /* qnew real + slot qnew for s */
  STEMP=qnew;
  for(i=0;i<qold;i++){
    int dst = (i==posS)? STEMP : o2n[i];
    if(dst<0) continue;
    if(omate[i]==-1) bmate[dst]=-1;           /* inner */
    else if(omate[i]>=0){                       /* outer: translate partner */
      int op=omate[i]; int pdst = (op==posS)? STEMP : o2n[op];
      bmate[dst]= pdst>=0? pdst : -2;
    } /* else bare stays -2 */
  }
}

int main(int argc,char*argv[]){
  if(argc<3){fprintf(stderr,"usage: dual m n [v]\n");return 1;}
  m=atoi(argv[1]); n=atoi(argv[2]); int verbose=argc>3; build();
  int q=frontier_before(0),i; for(i=0;i<q;i++) mate[i]=-2;
  unsigned char key[MAXF]; int kl=keyof(q,key);
  curkp=malloc(1<<20); curoff=malloc(sizeof(long)); curkl=malloc(sizeof(int)); curw=malloc(sizeof(u64));
  memcpy(curkp,key,kl); curoff[0]=0; curkl[0]=kl; curw[0]=1; ncur=1;

  int s;
  for(s=0;s<V;s++){
    int qold=frontier_before(s); posS=ifrb[s];
    static int frold[MAXF]; for(i=0;i<qold;i++) frold[i]=fr[i];
    qnew=frontier_before(s+1); static int ifrnew[64*64+2], frnew[MAXF];
    for(i=0;i<qnew;i++) frnew[i]=fr[i]; for(i=0;i<=V;i++) ifrnew[i]=-1; for(i=0;i<qnew;i++) ifrnew[frnew[i]]=i;
    apexnew=ifrnew[V];
    int nbr[16],rr=0; nbr[rr++]=apexnew;
    for(i=0;i<ND[s];i++){ int w=NB[s][i]; if(w>s && ifrnew[w]>=0) nbr[rr++]=ifrnew[w]; }
    for(i=0;i<qold;i++) o2n[i]= (frold[i]==s)? -1 : ifrnew[frold[i]];

    nr=0; kuse=0;
    long si;
    for(si=0;si<ncur;si++){
      unsigned char*ok=curkp+curoff[si]; int okl=curkl[si]; u64 w=curw[si];
      static int omate[MAXF];
      for(i=0;i<okl;i++){int c=ok[i]; omate[i]= c==0?-2 : c==255?-1 : c-1;}
      int deg = omate[posS]==-2? 0 : omate[posS]==-1? 2 : 1;   /* s's current degree */
      int need = 2-deg;
      unsigned char nk[MAXF]; int nl,a,b;
      if(need==0){
        build_bmate(omate,qold); for(i=0;i<qnew;i++) mate[i]=bmate[i];
        nl=keyof(qnew,nk); emit(nk,nl,w);
      } else if(need==1){
        for(a=0;a<rr;a++){
          build_bmate(omate,qold); for(i=0;i<=qnew;i++) mate[i]=bmate[i];
          if(add_derived(STEMP,nbr[a])){ nl=keyof(qnew,nk); emit(nk,nl,w); }
          else if(cycle){ int mp=completion_mp(qnew,apexnew,frnew,s); if(mp) cnt[mp]+=w; }
        }
      } else { /* need==2 */
        for(a=0;a<rr;a++)for(b=a+1;b<rr;b++){
          build_bmate(omate,qold); for(i=0;i<=qnew;i++) mate[i]=bmate[i];
          if(!add_derived(STEMP,nbr[a])){ if(cycle){ int mp=completion_mp(qnew,apexnew,frnew,s); if(mp) cnt[mp]+=w; } continue; }
          if(add_derived(STEMP,nbr[b])){ nl=keyof(qnew,nk); emit(nk,nl,w); }
          else if(cycle){ int mp=completion_mp(qnew,apexnew,frnew,s); if(mp) cnt[mp]+=w; }
        }
      }
    }
    cb=kp; qsort(rc,nr,sizeof(Rec),rcmp);
    curoff=realloc(curoff,(nr+1)*sizeof(long)); curkl=realloc(curkl,(nr+1)*sizeof(int));
    curw=realloc(curw,(nr+1)*sizeof(u64)); curkp=realloc(curkp,kuse+1);
    long k=0,use=0;
    for(si=0;si<nr;){ long j=si+1; u64 sw=rc[si].w;
      while(j<nr&&rc[j].len==rc[si].len&&memcmp(kp+rc[j].off,kp+rc[si].off,rc[si].len)==0){sw+=rc[j].w;j++;}
      memcpy(curkp+use,kp+rc[si].off,rc[si].len); curoff[k]=use; curkl[k]=rc[si].len; curw[k]=sw; use+=rc[si].len; k++; si=j; }
    ncur=k;
    if(verbose) fprintf(stderr,"s=%d col%d r%d states=%ld\n",s,s/m,s%m,ncur);
  }
  for(i=m;i<=V;i+=m) if(cnt[i]) printf("open %dx%d = %llu\n",m,i/m,cnt[i]);
  return 0;
}
