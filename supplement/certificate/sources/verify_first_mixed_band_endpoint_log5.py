"""Rigorous Arb--Feshbach certificate at a=log(5)/2.

Certified parameter sets:
    python verify_first_mixed_band_endpoint_log5.py 256 64 160 +
    python verify_first_mixed_band_endpoint_log5.py 256 60 120 -
"""
from fractions import Fraction
from flint import arb, acb, arb_mat, ctx
import sys,time
PAR=sys.argv[4] if len(sys.argv)>4 else '+'
ctx.prec=int(sys.argv[1]) if len(sys.argv)>1 else 256
R=int(sys.argv[2]) if len(sys.argv)>2 else (64 if PAR=='+' else 60)
N=int(sys.argv[3]) if len(sys.argv)>3 else (160 if PAR=='+' else 120)
if PAR not in {'+','-'}:
    raise ValueError('parity must be + or -')

def A(x):
    x=Fraction(x); return arb(x.numerator)/x.denominator
pi=arb.pi(); a=arb(5).log()/2
eps=1 if PAR=='+' else -1

def rate_value(k): return acb(A(k[0]),A(k[1])*pi/a)
def negkey(k): return (-k[0],-k[1])
def Iexp(k,lo,hi):
    if k==(Fraction(0),Fraction(0)): return acb(hi-lo)
    z=rate_value(k); return ((z*hi).exp()-(z*lo).exp())/z
def inner(r1,r2,lo,hi):
    s=acb(0)
    for c1,k1 in r1:
      for c2,k2 in r2:
        k=(k1[0]+k2[0],-k1[1]+k2[1]); s+=c1.conjugate()*c2*Iexp(k,lo,hi)
    return s.real
def basis(k):
    if PAR=='+':
      if k==0:return [(acb(1/a.sqrt()),(Fraction(0),Fraction(0)))]
      n=(arb(2)/a).sqrt();return [(acb(n/2),(Fraction(0),Fraction(k))),(acb(n/2),(Fraction(0),Fraction(-k)))]
    q=Fraction(2*k+1,2);n=(arb(2)/a).sqrt();return [(acb(0,-n/2),(Fraction(0),q)),(acb(0,n/2),(Fraction(0),-q))]
def transform(rep,kind,d):
    out=[]
    for c,k in rep:
      z=rate_value(k)
      if kind=='ref':out.append((c*(z*d).exp(),negkey(k)))
      elif kind=='plus':out.append((c*(z*d).exp(),k))
      elif kind=='minus':out.append((c*(-z*d).exp(),k))
    return out
def rankrep(n):
    ck=Fraction(4*n+1,2);c=A(ck);em=(-c*a).exp()
    if PAR=='+':
      den=(1-em*em).sqrt();return [(acb(em/den),(ck,Fraction(0))),(acb(em/den),(-ck,Fraction(0)))]
    den=(1+em*em).sqrt();return [(acb(em/den),(ck,Fraction(0))),(acb(-em/den),(-ck,Fraction(0)))]
def polar():return [(acb(1),(Fraction(1,2),Fraction(0))),(acb(1 if PAR=='+' else -1),(Fraction(-1,2),Fraction(0)))]
def scols(M,sg):return arb_mat([[M[i,j]*sg[j] for j in range(M.ncols())] for i in range(M.nrows())])

# Operator terms coeff*T, each T a partial reflection or translation.
d2=arb(2).log();d3=arb(3).log();d4=2*d2
alpha2=d2/arb(2).sqrt();alpha3=d3/arb(3).sqrt();alpha4=d2/2
TERMS=[]
# n=2, d2<a: -alpha(S+S*) - eps alpha J
TERMS += [(-eps*alpha2,'ref',d2,arb(0),d2),(-alpha2,'plus',d2,arb(0),a-d2),(-alpha2,'minus',d2,d2,a)]
# n=3,4: pure cross-centre reflection
TERMS += [(-eps*alpha3,'ref',d3,d3-a,a),(-eps*alpha4,'ref',d4,d4-a,a)]
# norm bound: ||S+S*||=1 since d2>a/2; each reflection norm=1
prime_bound=2*alpha2+alpha3+alpha4

start=time.time();BAS=[basis(i) for i in range(R)];FUN=[rankrep(n) for n in range(N)]+[polar()];M=N+1;sg=[1]*N+[1 if PAR=='+' else -1]
H=arb_mat([[inner(BAS[i],FUN[j],arb(0),a) for j in range(M)] for i in range(R)])
K=arb_mat([[inner(FUN[i],FUN[j],arb(0),a) for j in range(M)] for i in range(M)])
Prows=[[arb(0) for _ in range(R)] for _ in range(R)];Grows=[[arb(0) for _ in range(M)] for _ in range(R)]
TB=[]
for coeff,kind,d,lo,hi in TERMS:
  rb=[transform(x,kind,d) for x in BAS];rf=[transform(x,kind,d) for x in FUN];TB.append((coeff,lo,hi,rb))
  for i in range(R):
    for j in range(R):Prows[i][j]+=coeff*inner(BAS[i],rb[j],lo,hi)
    for j in range(M):Grows[i][j]+=coeff*inner(BAS[i],rf[j],lo,hi)
P=arb_mat(Prows);G=arb_mat(Grows)
P2r=[[arb(0) for _ in range(R)] for _ in range(R)]
for c1,l1,h1,r1 in TB:
 for c2,l2,h2,r2 in TB:
  lo=l1 if l1>l2 else l2;hi=h1 if h1<h2 else h2
  if not (hi>lo):continue
  for i in range(R):
   for j in range(R):P2r[i][j]+=c1*c2*inner(r1[i],r2[j],lo,hi)
P2=arb_mat(P2r);HS=scols(H,sg);GS=scols(G,sg)
Bpp=P+HS*H.transpose();B2=P2+GS*H.transpose()+HS*G.transpose()+HS*K*HS.transpose();CC=B2-Bpp*Bpp
D=[]
for k in range(R):
 q=Fraction(k) if PAR=='+' else Fraction(2*k+1,2);kap=A(q)*pi/a;D.append(acb(arb(1)/4,kap/2).digamma().real-pi.log())
q=Fraction(R) if PAR=='+' else Fraction(2*R+1,2);kap=A(q)*pi/a;dr=acb(arb(1)/4,kap/2).digamma().real-pi.log();mu=dr-prime_bound
if PAR=='-':mu-=2*(a.sinh()-a)
print('mu',mu)
S=arb_mat([[Bpp[i,j]-CC[i,j]/mu+(D[i] if i==j else 0) for j in range(R)] for i in range(R)]);S=(S+S.transpose())/2
L=[[arb(0) for _ in range(R)] for _ in range(R)];dd=[arb(0) for _ in range(R)];mn=None
for i in range(R):
 p=S[i,i]
 for k in range(i):p-=L[i][k]*L[i][k]*dd[k]
 if not (p>0):print('FAIL',i,p);sys.exit(2)
 dd[i]=p;lo=p.lower();mn=lo if mn is None or lo<mn else mn;L[i][i]=arb(1)
 for j in range(i+1,R):
  z=S[j,i]
  for k in range(i):z-=L[j][k]*L[i][k]*dd[k]
  L[j][i]=z/p
print('CERTIFIED',PAR,'minpivot',mn,'time',time.time()-start)
