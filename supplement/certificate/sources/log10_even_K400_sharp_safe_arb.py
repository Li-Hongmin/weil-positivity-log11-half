"""Safe sharp-weight K400 joint-covariance gate for either parity."""
from fractions import Fraction
import numpy as np,time
from flint import arb,acb,arb_mat,ctx
ARGV=__import__('sys').argv
FULL400='--full400' in ARGV;CORE7='--core7' in ARGV;CLEAN_NO_YUKAWA='--clean-no-yukawa' in ARGV;LOG11='--log11' in ARGV;ODD='--odd' in ARGV
PAR='-' if ODD else '+';PTAG='odd' if ODD else 'even';EPS=-1 if ODD else 1;FULL_CORE=6 if ODD else 7;CORE=7 if CORE7 else 1
if ODD and not (FULL400 and CLEAN_NO_YUKAWA):raise ValueError('odd gate currently requires --full400 --clean-no-yukawa')
FAMILY_PATH=None;FAMILY_M=0
if '--family' in ARGV:
 fi=ARGV.index('--family');FAMILY_PATH=ARGV[fi+1];FAMILY_M=int(ARGV[fi+2])
TAG='log11' if LOG11 else 'log10';ANUM=(11 if LOG11 else 10);ADEN=1;PREC=192
if '--cache-prefix' in ARGV:
 ci=ARGV.index('--cache-prefix');TAG=ARGV[ci+1];PREC=int(ARGV[ci+2]);ANUM=int(ARGV[ci+3]);ADEN=int(ARGV[ci+4])
INPREC=PREC
if '--input-precision' in ARGV:
 ii=ARGV.index('--input-precision');INPREC=int(ARGV[ii+1])
QF=2000
if '--qf' in ARGV:
 qi=ARGV.index('--qf');QF=int(ARGV[qi+1])
if QF<2000:raise ValueError('--qf must be at least 2000')
PRIME_BOUND_OVERRIDE=None
if '--prime-bound' in ARGV:
 pi=ARGV.index('--prime-bound');PRIME_BOUND_OVERRIDE=ARGV[pi+1]
FREEZE_WITNESS=None
if '--freeze-negative-witness' in ARGV:
 wi=ARGV.index('--freeze-negative-witness')
 if wi+1>=len(ARGV):raise ValueError('--freeze-negative-witness requires a path')
 FREEZE_WITNESS=ARGV[wi+1]
 if not (FULL400 and CLEAN_NO_YUKAWA and FAMILY_PATH is not None):
  raise ValueError('freezing requires --full400 --clean-no-yukawa --family')
TAIL_PROBE_END=None
if '--tail-probe-end' in ARGV:
 ti=ARGV.index('--tail-probe-end')
 if ti+1>=len(ARGV):raise ValueError('--tail-probe-end requires an integer')
 TAIL_PROBE_END=int(ARGV[ti+1])
 if FREEZE_WITNESS is None:raise ValueError('tail probe requires --freeze-negative-witness')
 if TAIL_PROBE_END<=QF:raise ValueError('tail probe endpoint must exceed --qf')
ctx.prec=PREC;AVAL=arb(ANUM)/arb(ADEN);E=400;R=280;M=CORE+120;N=500;Q0=400;CH=100;a=AVAL.log()/2;pi=arb.pi();prime_bound=arb(PRIME_BOUND_OVERRIDE) if PRIME_BOUND_OVERRIDE is not None else arb('2.978071' if TAG=='log10' else '2.931158');rho=prime_bound+(2*(a.sinh()-a) if ODD else arb(0));t0=time.time()
if FULL400:M=400
def A(x):x=Fraction(x);return arb(x.numerator)/x.denominator
def rate(k):return acb(A(k[0]),A(k[1])*pi/a)
def Iexp(k,lo,hi):
 if k==(Fraction(0),Fraction(0)):return acb(hi-lo)
 z=rate(k);return ((z*hi).exp()-(z*lo).exp())/z
def inner(r1,r2,lo,hi):
 s=acb(0)
 for c1,k1 in r1:
  for c2,k2 in r2:s+=c1.conjugate()*c2*Iexp((k1[0]+k2[0],-k1[1]+k2[1]),lo,hi)
 return s.real
def basis(k):
 if not ODD:
  if k==0:return [(acb(1/a.sqrt()),(Fraction(0),Fraction(0)))]
  c=(arb(2)/a).sqrt()/2;return [(acb(c),(Fraction(0),Fraction(k))),(acb(c),(Fraction(0),Fraction(-k)))]
 q=Fraction(2*k+1,2);c=(arb(2)/a).sqrt()/2
 return [(acb(0,-c),(Fraction(0),q)),(acb(0,c),(Fraction(0),-q))]
def transform(rep,kind,d):
 out=[]
 for c,k in rep:
  z=rate(k)
  if kind=='ref':out.append((c*(z*d).exp(),(-k[0],-k[1])))
  elif kind=='plus':out.append((c*(z*d).exp(),k))
  else:out.append((c*(-z*d).exp(),k))
 return out
def rankrep(n):
 ck=Fraction(4*n+1,2);c=A(ck);em=(-c*a).exp();den=(1+em*em).sqrt() if ODD else (1-em*em).sqrt()
 return [(acb(em/den),(ck,Fraction(0))),(acb((-em if ODD else em)/den),(-ck,Fraction(0)))]
def polar():return [(acb(1),(Fraction(1,2),Fraction(0))),(acb(-1 if ODD else 1),(Fraction(-1,2),Fraction(0)))]
def rank_overlap(rep,n):
 c=arb(2*n)+arb('0.5');em=(-c*a).exp();den=(1+em*em).sqrt() if ODD else (1-em*em).sqrt();ans=acb(0)
 for cb,k in rep:
  w=A(k[1])*pi/a;iw=acb(0,w)
  right=((acb(0,-w*a)).exp()-em)/(c-iw);left=em*(1-(-(c+iw)*a).exp())/(c+iw)
  ans+=cb.conjugate()*(right-left if ODD else right+left)/den
 return ans.real
def fover(rep,n):return rank_overlap(rep,n) if n<N else inner(rep,polar(),arb(0),a)
BE=[basis(k) for k in range(E)];FUN=[rankrep(n) for n in range(N)]+[polar()];HE=arb_mat([[fover(BE[i],n)*(EPS if n==N else 1) for n in range(N+1)] for i in range(E)])
d2=arb(2).log();d3=arb(3).log();d4=2*d2;d5=arb(5).log();d7=arb(7).log();d8=3*d2;d9=2*d3;d11=arb(11).log()
aa={2:d2/arb(2).sqrt(),3:d3/arb(3).sqrt(),4:d2/2,5:d5/arb(5).sqrt(),7:d7/arb(7).sqrt(),8:d2/arb(8).sqrt(),9:d3/3,11:d11/arb(11).sqrt()}
terms=[(-EPS*aa[2],'ref',d2,arb(0),d2),(-aa[2],'plus',d2,arb(0),a-d2),(-aa[2],'minus',d2,d2,a),(-EPS*aa[3],'ref',d3,arb(0),d3),(-aa[3],'plus',d3,arb(0),a-d3),(-aa[3],'minus',d3,d3,a)]
for n,d in [(4,d4),(5,d5),(7,d7),(8,d8),(9,d9),(11,d11)]:terms.append((-EPS*aa[n],'ref',d,d-a,a))
# Projection Z: first column fixed normalized physical f; remaining Q1 identity.
if not FULL400:
 c=np.load('scratch/log10_even_post_yukawa_physical_graph_arb_256.npz');fv=[arb(str(v)) for v in c['f']]
 ZC=arb_mat([[arb(str(c['Znew'][i,j])) for j in range(7)] for i in range(R)]) if CORE7 else arb_mat([[fv[i]] for i in range(R)])
def project(Mat):
 if FULL400:return Mat
 out=arb_mat(M,M)
 PP=arb_mat([[Mat[i,j] for j in range(R)] for i in range(R)]);PQ=arb_mat([[Mat[i,R+j] for j in range(120)] for i in range(R)]);cc=ZC.transpose()*PP*ZC;cq=ZC.transpose()*PQ
 for i in range(CORE):
  for j in range(CORE):out[i,j]=cc[i,j]
  for j in range(120):out[i,CORE+j]=cq[i,j];out[CORE+j,i]=cq[i,j]
 for j in range(120):
  for k in range(120):out[CORE+j,CORE+k]=Mat[R+j,R+k]
 return out
BASE_PATH=(f'scratch/{TAG}_{PTAG}_K400_N500_clean_full_contract_base_arb_{INPREC}.npz' if CLEAN_NO_YUKAWA else 'scratch/log10_even_K400_full_contract_base_arb_192.npz') if FULL400 else ('scratch/log10_even_K400_projected_core7_arb_192.npz' if CORE7 else 'scratch/log10_even_K400_projected_base_arb_192.npz')
base=np.load(BASE_PATH,allow_pickle=FULL400)
if FULL400 and CLEAN_NO_YUKAWA and 'parity' in base.files:
 assert str(base['parity'])==PAR and str(base['cache_prefix'])==TAG and int(base['precision'])==INPREC
 assert int(base['N'])==N and int(base['R'])==R and bool(base['clean'])
 assert (arb(str(base['a']))-a).contains(0)
def load400(raw):
 out=arb_mat(E,E)
 for i0,i1,j0,j1 in ((0,R,0,R),(0,R,R,E),(R,E,R,E)):
  b=arb_mat([[arb(str(raw[i,j])) for j in range(j0,j1)] for i in range(i0,i1)])
  for i in range(i1-i0):
   for j in range(j1-j0):out[i0+i,j0+j]=b[i,j];out[j0+j,i0+i]=b[i,j]
 return out
P0=load400(base['S']) if FULL400 else arb_mat([[arb(str(base['S'][i,j])) for j in range(M)] for i in range(M)])
RP=load400(base['Rcov']) if FULL400 else arb_mat([[arb(str(base['Rcov'][i,j])) for j in range(M)] for i in range(M)])
G1=arb_mat(M,M);R1=arb_mat(M,M)
GRAPH_PREC=INPREC if TAG!='log10' or ODD else 256
GRAPH_SOURCE=(f'{TAG}_{PTAG}_R280_N500_verified_graph_arb_{GRAPH_PREC}.npz'
              if CLEAN_NO_YUKAWA else 'log10_even_R280_N400_exact_small_map_256.npz')
ex=np.load('scratch/'+GRAPH_SOURCE) if FULL400 else None
whist=(ex['whist'] if CLEAN_NO_YUKAWA else ex['whist_arb']) if FULL400 else None
if FULL400 and CLEAN_NO_YUKAWA:
 assert int(ex['N'])==N and int(ex['R'])==R and str(ex['parity'])==PAR and str(ex['cache_prefix'])==TAG
 assert (arb(str(ex['a']))-a).contains(0)
def left_contract(W):
 if not FULL400:return W
 for row in whist:
  h=arb_mat([[arb(str(row[i]))] for i in range(R)]);hh=(h.transpose()*h)[0,0]
  top=arb_mat([[W[i,j] for j in range(W.ncols())] for i in range(R)])
  top-=h*((arb(2)/hh)*(h.transpose()*top))
  for i in range(R):
   for j in range(W.ncols()):W[i,j]=top[i,j]
 return W
if FAMILY_PATH is not None:
 if not (FULL400 and CLEAN_NO_YUKAWA):raise ValueError('family joint gate requires full400 clean base')
 fam=np.load(FAMILY_PATH);assert (str(fam['parity'])==PAR and int(fam['R'])==E and int(fam['N'])==N and int(fam['atom_lo'])==500)
 assert str(fam['cache_prefix'])==TAG and int(fam['precision'])==INPREC and bool(fam['exact_parseval'])
 assert 1<=FAMILY_M<=int(fam['atom_hi'])-int(fam['atom_lo'])
 # Materialise each compressed NPZ member once.  Re-indexing ``fam['b']``
 # inside the nested loop makes NumPy decompress the whole member on every
 # entry; caching the arrays changes only I/O complexity, not any Arb value.
 fam_b=fam['b'];fam_cross=fam['cross'];fam_gram_q=fam['gram_q']
 BF=arb_mat([[arb(str(fam_b[i,j])) for j in range(FAMILY_M)] for i in range(E)])
 XF=arb_mat([[arb(str(fam_cross[i,j])) for j in range(FAMILY_M)] for i in range(E)])
 GQ=arb_mat([[arb(str(fam_gram_q[i,j])) for j in range(FAMILY_M)] for i in range(FAMILY_M)])
 BF=left_contract(BF);XF=left_contract(XF)
 RTOT=RP+XF*BF.transpose()+BF*XF.transpose()+BF*GQ*BF.transpose()
 print('JOINT_FULL_YUKAWA_FAMILY',500,500+FAMILY_M,flush=True)
def tq_at(k):
 q=Fraction(k) if not ODD else Fraction(2*k+1,2)
 return acb(arb(1)/4,A(q)*pi/(2*a)).digamma().real-pi.log()-rho
if FAMILY_PATH is not None:
 assert (arb(str(fam['a']))-a).contains(0)
 assert (arb(str(fam['t400']))-tq_at(400)).contains(0)
for lo in range(Q0,QF,CH):
 hi=min(lo+CH,QF);BQ=[basis(q) for q in range(lo,hi)];HQ=arb_mat([[fover(BQ[j-lo],n) for n in range(N+1)] for j in range(lo,hi)])
 W=arb_mat(E,hi-lo)
 for co,kind,d,xlo,xhi in terms:
  TB=[transform(x,kind,d) for x in BQ]
  for i in range(E):
   for j in range(hi-lo):W[i,j]+=co*inner(BE[i],TB[j],xlo,xhi)
 W+=HE*HQ.transpose()
 if FULL400:
  U=left_contract(W)
  if FAMILY_PATH is not None:
   VQ=arb_mat([[rank_overlap(BQ[j],500+n) for j in range(hi-lo)] for n in range(FAMILY_M)])
   U+=BF*VQ
  Ws=arb_mat([[U[i,j]/tq_at(lo+j).sqrt() for j in range(hi-lo)] for i in range(M)]);R1+=U*U.transpose();G1+=Ws*Ws.transpose()
 else:
  U=arb_mat(M,hi-lo);UC=ZC.transpose()*arb_mat([[W[i,j] for j in range(hi-lo)] for i in range(R)])
  for j in range(hi-lo):
   for i in range(CORE):U[i,j]=UC[i,j]
   for i in range(120):U[CORE+i,j]=W[R+i,j]
  for j,q in enumerate(range(lo,hi)):
   tq=tq_at(q)
   for i in range(M):
    for k in range(i,M):
     vv=U[i,j]*U[k,j];R1[i,k]+=vv;G1[i,k]+=vv/tq
  for i in range(M):
   for k in range(i):R1[i,k]=R1[k,i];G1[i,k]=G1[k,i]
 print('ARB_SHARP_PREFIX_END',hi,'SECONDS',time.time()-t0,flush=True)
t400=tq_at(400);tQF=tq_at(QF)
if FAMILY_PATH is not None:
 G=P0+RP/t400+BF*BF.transpose()-G1-(RTOT-R1)/tQF
else:
 G=P0+RP*(1/t400-1/tQF)-G1+R1/tQF
if CLEAN_NO_YUKAWA:
 print('CLEAN_NO_ISOLATED_YUKAWA',True,flush=True)
else:
 # Add compressed rational PSD Yukawa T0 through its 121 projection.
 # This historical branch is not operator-promotable; the clean gate skips it.
 nu=arb(N)+arb('0.25');al=[arb(k)*pi/(2*a) for k in range(E)];F=[arb(0)]*E;Fp=[arb(0)]*E;F[0]=acb(nu).polygamma(1).real
 for k in range(1,E):
  zz=acb(nu,al[k]);p=zz.digamma();tr=zz.polygamma(1);F[k]=p.imag/al[k];Fp[k]=(al[k]*tr.real-p.imag)/(al[k]*al[k])
 D=[[arb(0) for _ in range(E)] for _ in range(E)];D[0][0]=F[0]/4
 for i in range(E):
  if i>0:D[i][i]=(F[i]+al[i]*Fp[i]/2)/4
  for j in range(i+1,E):
   v=F[j]/4 if i==0 else (al[j]*al[j]*F[j]-al[i]*al[i]*F[i])/(4*(al[j]*al[j]-al[i]*al[i]));D[i][j]=v;D[j][i]=v
 fac=arb(1)-(-nu*arb(100).log()).exp();rt2=arb(2).sqrt();T=arb_mat([[(arb(2)/a)*fac*((-1) if (i+j)%2 else 1)*D[i][j]/((rt2 if i==0 else arb(1))*(rt2 if j==0 else arb(1))) for j in range(E)] for i in range(E)])
 if FULL400:
  # T is still in physical coordinates; contract it once, blockwise.
  for row in whist:
   h=arb_mat([[arb(str(row[i]))] if i<R else [arb(0)] for i in range(E)]);hh=(h.transpose()*h)[0,0];Th=T*h
   T=T-h*((arb(2)/hh)*(h.transpose()*T))-Th*((arb(2)/hh)*h.transpose())+h*((arb(4)*(h.transpose()*Th)[0,0]/(hh*hh))*h.transpose())
  G+=T
 else:G+=project(T)
if FULL400:
 CORE=FULL_CORE;H=arb_mat([[G[i+CORE,j+CORE] for j in range(E-CORE)] for i in range(E-CORE)]);Y=arb_mat([[G[i+CORE,j] for j in range(CORE)] for i in range(E-CORE)])
else:H=arb_mat([[G[i+CORE,j+CORE] for j in range(120)] for i in range(120)]);Y=arb_mat([[G[i+CORE,j] for j in range(CORE)] for i in range(120)])
# Strict LDL certifies H>0 before verified solve.
NH=E-CORE if FULL400 else 120;L=[[arb(0) for _ in range(NH)] for _ in range(NH)];minp=None
for i in range(NH):
 p=H[i,i]-sum((L[i][k]*L[i][k] for k in range(i)),arb(0))
 if not (p>0):raise ArithmeticError(('H_cholesky',i,p))
 L[i][i]=p.sqrt();minp=p if minp is None or p.lower()<minp.lower() else minp
 for j in range(i+1,NH):L[j][i]=(H[j,i]-sum((L[j][k]*L[i][k] for k in range(i)),arb(0)))/L[i][i]
sol=H.solve(Y);SC=arb_mat([[G[i,j] for j in range(CORE)] for i in range(CORE)])-Y.transpose()*sol;SC=(SC+SC.transpose())/2
print('H_CHOLESKY_MIN_PIVOT',minp)
if not CORE7 and not FULL400:
 ray=SC[0,0];print('K400_SHARP_SAFE_DIRECTION_RAY',ray);print('RAY_OVER_1E12',ray/arb('1e-12'));np.savez('scratch/log10_even_K400_sharp_safe_direction_arb_192.npz',ray=np.array(str(ray)),Hmin=np.array(str(minp)))
else:
 Tld=[[SC[i,j] for j in range(CORE)] for i in range(CORE)];perm=list(range(CORE));pivs=[];negative_step=None
 for k in range(CORE):
  p=max(range(k,CORE),key=lambda i:float(Tld[i][i].mid()))
  if p!=k:
   Tld[k],Tld[p]=Tld[p],Tld[k]
   for row in Tld:row[k],row[p]=row[p],row[k]
  perm[k],perm[p]=perm[p],perm[k]
  d=Tld[k][k];pivs.append(d);print('CORE_PIVOT',k,d)
  if d.upper()<0:
   negative_step=k;break
  if not d.lower()>0:break
  for i in range(k+1,CORE):
   for j in range(i,CORE):
    v=Tld[j][i]-Tld[j][k]*Tld[i][k]/d;Tld[j][i]=v;Tld[i][j]=v
 print('CORE_PIVOT_PERM',perm,'PIVOTS',len(pivs))
 if FREEZE_WITNESS is not None:
  if negative_step is None:raise ArithmeticError(('no_strict_negative_core_pivot',pivs,perm))
  # This is the post-ordered-Householder coordinate system of G.  Keep the
  # rounded vector in that system: the tail U_q in this gate is defined there.
  inds=perm[:negative_step];last=perm[negative_step]
  Ap=arb_mat([[SC[i,j] for j in inds] for i in inds]);bp=arb_mat([[SC[i,last]] for i in inds])
  zp=Ap.solve(bp);z=arb_mat(CORE,1);z[last,0]=arb(1)
  for j,i in enumerate(inds):z[i,0]=-zp[j,0]
  xraw=arb_mat(E,1)
  for i in range(CORE):xraw[i,0]=z[i,0]
  yz=sol*z
  for i in range(E-CORE):xraw[CORE+i,0]=-yz[i,0]
  xunit=xraw/(xraw.transpose()*xraw)[0,0].sqrt()
  # ``mid_rad_10exp`` exposes a decimal midpoint as integer times a power of
  # ten.  The conversion below rounds that rational using only integer
  # arithmetic; no binary float chooses a certified witness coordinate.
  def nearest_dyadic(v,bits=224):
   mant,unused_rad,tenexp=v.mid_rad_10exp(100)
   mant=int(mant);tenexp=int(tenexp)
   if tenexp>=0:return mant*(10**tenexp)*(2**bits)
   den10=10**(-tenexp);num=mant*(2**bits)
   sign=1 if num>=0 else -1;num=abs(num)
   return sign*((2*num+den10)//(2*den10))
  WITNESS_BITS=224;den=2**WITNESS_BITS;nums=[nearest_dyadic(xunit[i,0],WITNESS_BITS) for i in range(E)]
  x=arb_mat([[arb(v)/den] for v in nums]);norm2=(x.transpose()*x)[0,0]
  quadratic=(x.transpose()*G*x)[0,0];ray=quadratic/norm2
  print('FROZEN_NEGATIVE_WITNESS_NORM2',norm2);print('FROZEN_NEGATIVE_WITNESS_RAY',ray)
  if not ray.upper()<0:raise ArithmeticError(('dyadic_full_witness_not_negative',ray))
  import hashlib
  def sha256(path):
   h=hashlib.sha256()
   with open(path,'rb') as fh:
    for block in iter(lambda:fh.read(1<<20),b''):h.update(block)
   return h.hexdigest()
  d_prefix=arb(0)
  if TAIL_PROBE_END is not None:
   # A direct positive prefix of D_QF in the *same* contracted coordinates.
   # It is a lower bound for D_QF(x), so if it exceeds -ray this particular
   # lower-gate witness is already repaired by the true Schur tail.
   for lo in range(QF,TAIL_PROBE_END,CH):
    hi=min(lo+CH,TAIL_PROBE_END);BQ=[basis(q) for q in range(lo,hi)];HQ=arb_mat([[fover(BQ[j-lo],n) for n in range(N+1)] for j in range(lo,hi)])
    W=arb_mat(E,hi-lo)
    for co,kind,d,xlo,xhi in terms:
     TB=[transform(v,kind,d) for v in BQ]
     for i in range(E):
      for j in range(hi-lo):W[i,j]+=co*inner(BE[i],TB[j],xlo,xhi)
    W+=HE*HQ.transpose();U=left_contract(W)
    VQ=arb_mat([[rank_overlap(BQ[j],500+n) for j in range(hi-lo)] for n in range(FAMILY_M)]);U+=BF*VQ
    for j,q in enumerate(range(lo,hi)):
     cq=(x.transpose()*arb_mat([[U[i,j]] for i in range(E)]))[0,0]
     d_prefix+=(1/tQF-1/tq_at(q))*cq*cq
    print('FROZEN_NEGATIVE_WITNESS_D_PREFIX_END',hi,d_prefix/norm2,'OVER_DEFICIT',(d_prefix/norm2)/(-ray),flush=True)
  np.savez(FREEZE_WITNESS,numerator=np.array([str(v) for v in nums]),denominator=np.array(str(den)),bits=np.array(WITNESS_BITS),coordinate=np.array('post-ordered-Householder'),norm2=np.array(str(norm2)),quadratic=np.array(str(quadratic)),rayleigh=np.array(str(ray)),negative_step=np.array(negative_step),negative_pivot=np.array(str(pivs[negative_step])),perm=np.array(perm),qf=np.array(QF),t400=np.array(str(t400)),tqf=np.array(str(tQF)),a=np.array(str(a)),endpoint_value=np.array(str(AVAL)),parity=np.array(PAR),R=np.array(R),N=np.array(N),family_atom_lo=np.array(500),family_atom_hi=np.array(500+FAMILY_M),tail_probe_end=np.array(-1 if TAIL_PROBE_END is None else TAIL_PROBE_END),tail_probe_weighted_prefix=np.array(str(d_prefix/norm2)),base_path=np.array(BASE_PATH),family_path=np.array(FAMILY_PATH),householder_path=np.array('scratch/'+GRAPH_SOURCE),base_sha256=np.array(sha256(BASE_PATH)),family_sha256=np.array(sha256(FAMILY_PATH)),householder_sha256=np.array(sha256('scratch/'+GRAPH_SOURCE)),source_path=np.array(__file__),source_sha256=np.array(sha256(__file__)))
  print('FROZEN_NEGATIVE_WITNESS_SAVED',FREEZE_WITNESS)
 # A successful process is itself part of the certificate contract.  A
 # negative or interval-indeterminate pivot must never fall through to the
 # NPZ writer with exit status zero, and a truncated pivot list is a failure.
 if len(pivs) != CORE:
  raise ArithmeticError(('incomplete_core_pivots',len(pivs),CORE,pivs,perm))
 for i,d in enumerate(pivs):
  if not d.lower()>0:
   raise ArithmeticError(('nonpositive_or_indeterminate_core_pivot',i,d,pivs,perm))
 qtag='' if QF==2000 else f'_q{QF}';np.savez((f'scratch/{TAG}_{PTAG}_K400_N500_clean_sharp_family_500_{500+FAMILY_M}_full{E-CORE}{qtag}_arb_{PREC}.npz' if FAMILY_PATH is not None else f'scratch/{TAG}_{PTAG}_K400_N500_clean_sharp_full{E-CORE}{qtag}_arb_{PREC}.npz' if CLEAN_NO_YUKAWA else 'scratch/log10_even_K400_sharp_safe_full393_arb_192.npz') if FULL400 else 'scratch/log10_even_K400_sharp_safe_core7_arb_192.npz',S=np.array([[str(SC[i,j]) for j in range(CORE)] for i in range(CORE)]),Hmin=np.array(str(minp)),pivots=np.array([str(x) for x in pivs]),perm=np.array(perm),N=np.array(N),clean=np.array(CLEAN_NO_YUKAWA),parity=np.array(PAR),cache_prefix=np.array(TAG),endpoint_value=np.array(str(AVAL)),precision=np.array(PREC),qf=np.array(QF),t400=np.array(str(t400)),tqf=np.array(str(tQF)),base_source=np.array(str(base['source']) if 'source' in base.files else ''),family_source=np.array('' if FAMILY_PATH is None else FAMILY_PATH),family_m=np.array(FAMILY_M),householder_source=np.array(GRAPH_SOURCE))
