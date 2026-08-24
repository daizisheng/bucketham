#!/usr/bin/env python3
# CRT orchestration: run the SpMV under K primes, combine to exact integers.
import subprocess, sys, re
from functools import reduce

def is_prime(n):
    if n<2: return False
    for p in [2,3,5,7,11,13,17,19,23,29,31,37]:
        if n%p==0: return n==p
    d=n-1; r=0
    while d%2==0: d//=2; r+=1
    for a in [2,3,5,7,11,13,17,19,23,29,31,37]:
        x=pow(a,d,n)
        if x==1 or x==n-1: continue
        for _ in range(r-1):
            x=x*x%n
            if x==n-1: break
        else: return False
    return True

def primes_below(limit, k):
    out=[]; n=limit-1
    while len(out)<k:
        if is_prime(n): out.append(n)
        n-=1
    return out

def run(exe, tables, Nto, prime):
    r=subprocess.run([exe,'run',tables,str(Nto),str(prime)],capture_output=True,text=True)
    d={}
    for line in r.stdout.splitlines():
        mo=re.match(r'open (\d+)x(\d+) = (\d+)',line)
        if mo: d[int(mo.group(2))]=int(mo.group(3))
    return d

def crt(residues, primes):
    # residues[i] mod primes[i] -> x mod prod
    x=0; M=1
    for a,p in zip(residues,primes):
        # combine x (mod M) with a (mod p)
        g=pow(M,-1,p)
        x = x + M*((a-x)*g % p)
        M*=p
    return x%M

if __name__=='__main__':
    exe,tables,m,Nto,K = sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5])
    primes=primes_below(2**31, K)
    print(f"using {K} primes below 2^31, product ~2^{K*31}")
    tabs=[run(exe,tables,Nto,p) for p in primes]
    ref={}
    if m==5:
        try:
            for line in open('../data/S5_open.txt'):
                a=line.split(); ref[int(a[0])]=int(a[1])
        except: pass
    else:  # exact run (MODP=0) is a reference where the value fits in u64
        ex=run(exe,tables,Nto,0)
        for c,v in ex.items():
            if v < (1<<63): ref[c]=v   # trust only clearly-non-overflowed values
    cols=sorted(set().union(*[set(t) for t in tabs]))
    ok=True
    for c in cols:
        res=[t.get(c) for t in tabs]
        if any(r is None for r in res): continue
        x=crt(res,primes)
        tag=""
        if c in ref and x < (1<<64):   # exact run valid only when true value fits u64
            tag = " OK" if x==ref[c] else f" MISMATCH(ref={ref[c]})"
            if x!=ref[c]: ok=False
        elif x >= (1<<64):
            tag = "  (exact, >u64)"
        print(f"open {m}x{c} = {x}{tag}")
    print("ALL CRT VALUES MATCH REFERENCE" if ok else "SOME MISMATCH")
