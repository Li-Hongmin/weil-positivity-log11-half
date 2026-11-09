"""Executed inside the temporary log7 probe after its Arb matrix S is built."""

import sys
import numpy as np

assert (PAR == "+" and R == 128) or (R == 280 and N == 500)
BLOCK = 4
VECTOR_BITS = 44
CHOLESKY_BITS = 52
C0 = arb(1) / 2000000

B = S
Bmid = np.array([[float(B[i, j].mid()) for j in range(R)] for i in range(R)])
initial_eigs = np.linalg.eigvalsh(Bmid)[:8]
Ufloat = np.eye(R)
whist = []

# Four exact orthogonal Householder congruences.  At stage k the generator
# has its first k entries zero, so previously isolated coordinates stay fixed.
for k in range(BLOCK):
    trailing = (Bmid[k:, k:] + Bmid[k:, k:].T) / 2
    _, vecs = np.linalg.eigh(trailing)
    q = vecs[:, 0]
    sign = 1.0 if q[0] >= 0 else -1.0
    wf = np.zeros(R)
    wf[k:] = q
    wf[k] -= sign
    wint = np.rint(wf * 2.0**VECTOR_BITS).astype(object)
    w = [arb(int(z)) / (2**VECTOR_BITS) for z in wint]
    whist.append(w)
    t = sum((z*z for z in w), arb(0))
    assert t.lower() > 0
    sw = [sum((B[i, j]*w[j] for j in range(R)), arb(0)) for i in range(R)]
    wsw = sum((w[i]*sw[i] for i in range(R)), arb(0))
    two = arb(2)/t
    four = arb(4)*wsw/t**2
    B = arb_mat([[B[i, j]-two*(w[i]*sw[j]+sw[i]*w[j])+four*w[i]*w[j]
                  for j in range(R)] for i in range(R)])
    # Float update only selects the next witness.
    wfd = np.array([float(z.mid()) for z in w])
    Hf = np.eye(R)-2*np.outer(wfd,wfd)/np.dot(wfd,wfd)
    Bmid = Hf @ Bmid @ Hf
    Ufloat = Ufloat @ Hf

m = R-BLOCK
Cminus = arb_mat([[B[i+BLOCK,j+BLOCK]-(C0 if i==j else arb(0))
                   for j in range(m)] for i in range(m)])
Cmid = np.array([[float(Cminus[i,j].mid()) for j in range(m)] for i in range(m)])
try:
    Lf = np.linalg.cholesky((Cmid+Cmid.T)/2)
except np.linalg.LinAlgError:
    cache_prefix = globals().get('CACHE_PREFIX', 'log8')
    np.savez(
        f'scratch/{cache_prefix}_{"even" if PAR=="+" else "odd"}_R{R}_N{N}_tail_nonpd_arb_{ctx.prec}.npz',
        Cminus=np.array([[str(Cminus[i,j]) for j in range(m)] for i in range(m)]),
        Cmid=Cmid, R=np.array(R), N=np.array(N), BLOCK=np.array(BLOCK),
        C0=np.array(str(C0)), parity=np.array(PAR), precision=np.array(ctx.prec),
        cache_prefix=np.array(cache_prefix), a=np.array(str(a)), mu=np.array(str(mu)))
    print('EXPORTED_NONPD_TAIL', cache_prefix, R, N, ctx.prec, flush=True)
    raise
scale = 2**CHOLESKY_BITS
L = arb_mat([[arb(int(np.rint(Lf[i,j]*scale)))/scale if j<=i else arb(0)
              for j in range(m)] for i in range(m)])
res = Cminus-L*L.transpose()
resinf = arb(0)
for i in range(m):
    rowsum = sum((abs(res[i,j]) for j in range(m)),arb(0))
    if rowsum.upper()>resinf.upper(): resinf=rowsum
clower=C0-resinf
assert clower.lower()>0

A4=arb_mat([[B[i,j] for j in range(BLOCK)] for i in range(BLOCK)])
X=arb_mat([[B[i,j+BLOCK] for j in range(m)] for i in range(BLOCK)])
Cmat=arb_mat([[B[i+BLOCK,j+BLOCK] for j in range(m)] for i in range(m)])
# Once C>0 is independently certified, Arb's verified solve gives the sharp
# block Schur complement.  The scalar C-lower estimate is far too lossy at
# the 1e-16 near-null scale.
Y=Cmat.solve(X.transpose())
small=A4-X*Y
smallmid=np.array([[float(small[i,j].mid()) for j in range(BLOCK)] for i in range(BLOCK)])
# Export the exact five-dimensional block and its pullback columns for the
# endpoint selected by the wrapper.  Generated interval caches remain in
# scratch and are not theorem inputs until recreated by the documented run.
cache_prefix=globals().get('CACHE_PREFIX','log8')
export_case=((cache_prefix=='log8' and PAR=='-' and R==200)
             or (cache_prefix=='log9' and R==240))
if export_case and BLOCK==5 and N==400:
    zcols=[]
    for col in range(BLOCK):
        zz=[arb(1) if i==col else arb(0) for i in range(BLOCK)]
        zz += [-Y[j,col] for j in range(m)]
        for ww in reversed(whist):
            wt=sum((ww[i]*zz[i] for i in range(R)),arb(0));wn=sum((v*v for v in ww),arb(0))
            zz=[zz[i]-2*ww[i]*wt/wn for i in range(R)]
        zcols.append(zz)
    np.savez('scratch/'+cache_prefix+'_'+('even' if PAR=='+' else 'odd')+f'_R{R}_N400_exact_small_map.npz',
             small_arb=np.array([[str(small[i,j]) for j in range(BLOCK)] for i in range(BLOCK)]),
             z_arb=np.array([[str(zcols[j][i]) for j in range(BLOCK)] for i in range(R)]),
             x_arb=np.array([[str(X[i,j]) for j in range(m)] for i in range(BLOCK)]),
             c_arb=np.array([[str(Cmat[i,j]) for j in range(m)] for i in range(m)]),
             whist_arb=np.array([[str(w[i]) for i in range(R)] for w in whist]))
    if '--export-small-map-only' in sys.argv:
        print('EXPORTED_EXACT_SMALL_MAP',flush=True);raise SystemExit(0)
sew,sev=np.linalg.eigh((smallmid+smallmid.T)/2)
print('small_midpoint_eigs',sew,flush=True)
qint=np.rint(sev[:,0]*2**50).astype(object)
q=[arb(int(z))/2**50 for z in qint]
ray=sum((q[i]*small[i,j]*q[j] for i in range(BLOCK) for j in range(BLOCK)),arb(0))
qnorm=sum((z*z for z in q),arb(0))
print('small_negative_rayleigh',ray/qnorm,flush=True)

# Arb pullback of [q,-C^{-1}X^Tq] through the exact Householders.
ytail=[-sum((Y[j,i]*q[i] for i in range(BLOCK)),arb(0)) for j in range(m)]
xarb=q+ytail
for w in reversed(whist):
    wt=sum((w[i]*xarb[i] for i in range(R)),arb(0))
    ww=sum((z*z for z in w),arb(0))
    xarb=[xarb[i]-2*w[i]*wt/ww for i in range(R)]
xarb_norm=sum((z*z for z in xarb),arb(0)).sqrt()
xarb=[z/xarb_norm for z in xarb]
arb_sray=sum((xarb[i]*S[i,j]*xarb[j] for i in range(R) for j in range(R)),arb(0))
arb_ccray=sum((xarb[i]*CC[i,j]*xarb[j] for i in range(R) for j in range(R)),arb(0))
print('arb_pulled_schur_rayleigh',arb_sray,flush=True)
print('arb_pulled_CC_rayleigh',arb_ccray,flush=True)

# Export the exact Schur witness pulled back to the original Fourier core.
# This is scratch-only diagnostic data used by the directional-resolvent gate.
qf=np.array([float(z.mid()) for z in q])
Cfullmid=np.array([[float(Cmat[i,j].mid()) for j in range(m)] for i in range(m)])
Xfullmid=np.array([[float(X[i,j].mid()) for j in range(m)] for i in range(BLOCK)])
yfull=np.concatenate([qf,-np.linalg.solve(Cfullmid,Xfullmid.T@qf)])
xfull=np.array([float(z.mid()) for z in xarb])
np.savez('scratch/log8_'+('even' if PAR=='+' else 'odd')+
         f'_R{R}_N{N}_schur_witness.npz',x=xfull,
         x_arb=np.array([str(z.mid()) for z in xarb]),
         small=smallmid,small_eigs=sew,
         scalar_S=np.array([[float(S[i,j].mid()) for j in range(R)] for i in range(R)]),
         CC=np.array([[float(CC[i,j].mid()) for j in range(R)] for i in range(R)]),
         mu=float(mu.mid()))
print('exported_schur_witness_rayleigh',float(xfull@np.array(
      [[float(S[i,j].mid()) for j in range(R)] for i in range(R)])@xfull),flush=True)

# Fixed log8/R=200 odd diagnostic: symmetric diagonal-pivoted interval LDL.
# Explicit Schur updates avoid the severe dependency growth of the original
# unpivoted recurrence when a ~1e-27 pivot is encountered before large ones.
if PAR=='-' and R==200 and BLOCK==5:
    T=[[small[i,j] for j in range(BLOCK)] for i in range(BLOCK)]
    perm=list(range(BLOCK)); pp=[]
    for k in range(BLOCK):
        p=max(range(k,BLOCK),key=lambda i:float(T[i][i].mid()))
        if p!=k:
            T[k],T[p]=T[p],T[k]
            for row in T: row[k],row[p]=row[p],row[k]
            perm[k],perm[p]=perm[p],perm[k]
        pivk=T[k][k]
        pp.append(pivk)
        print('pivoted_small_pivot',k,pivk,flush=True)
        if not pivk.lower()>0:
            # Exact witness for the failed final Schur pivot, pulled back first
            # through the large C block and then through all Householders.
            assert k==BLOCK-1
            inds=perm[:k]; last=perm[k]
            Aneg=arb_mat([[small[i,j] for j in inds] for i in inds])
            bneg=arb_mat([[small[i,last]] for i in inds])
            yneg=Aneg.solve(bneg)
            qneg=[arb(0) for _ in range(BLOCK)];qneg[last]=arb(1)
            for ii,i0 in enumerate(inds):qneg[i0]=-yneg[ii,0]
            tailneg=[-sum((Y[j,i]*qneg[i] for i in range(BLOCK)),arb(0)) for j in range(m)]
            xneg=qneg+tailneg
            for ww in reversed(whist):
                wt=sum((ww[i]*xneg[i] for i in range(R)),arb(0));wn=sum((z*z for z in ww),arb(0))
                xneg=[xneg[i]-2*ww[i]*wt/wn for i in range(R)]
            xn=sum((z*z for z in xneg),arb(0)).sqrt();xneg=[z/xn for z in xneg]
            rneg=sum((xneg[i]*S[i,j]*xneg[j] for i in range(R) for j in range(R)),arb(0))
            print('pivot_failure_witness_rayleigh',rneg,flush=True)
            np.savez(f'scratch/log8_odd_R{R}_N{N}_pivot_failure_witness.npz',
                     x=np.array([float(z.mid()) for z in xneg]),
                     x_arb=np.array([str(z.mid()) for z in xneg]),rayleigh=str(rneg))
            print('PIVOTED_BLOCK_FAIL',k,pivk,'perm',perm,flush=True)
            break
        for i in range(k+1,BLOCK):
            for j in range(i,BLOCK):
                v=T[j][i]-T[j][k]*T[i][k]/pivk
                T[j][i]=v;T[i][j]=v
    else:
        print('PIVOTED_BLOCK_CERTIFIED','perm',perm,'minlower',
              min(z.lower() for z in pp),flush=True)
        raise SystemExit(0)

low=[[arb(0) for _ in range(BLOCK)] for _ in range(BLOCK)]
piv=[]
for i in range(BLOCK):
    p=small[i,i]
    for k in range(i): p-=low[i][k]**2*piv[k]
    print('small_pivot',i,p,flush=True)
    if not p.lower()>0:
        np.savez('scratch/log8_'+('even' if PAR=='+' else 'odd')+
                 f'_R{R}_N{N}_block_status.npz',certified=False,index=i,
                 pivots=np.array([str(z) for z in piv+[p]]))
        print('BLOCK_FAIL',BLOCK,i,p,flush=True)
        sys.exit(2)
    piv.append(p);low[i][i]=arb(1)
    for j in range(i+1,BLOCK):
        z=small[j,i]
        for k in range(i):z-=low[j][k]*low[i][k]*piv[k]
        low[j][i]=z/p

np.savez('scratch/log8_'+('even' if PAR=='+' else 'odd')+
         f'_R{R}_N{N}_block_status.npz',certified=True,index=BLOCK,
         pivots=np.array([str(z) for z in piv]))

print('initial_midpoint_eigs',initial_eigs,flush=True)
print('C0',C0,flush=True)
print('cholesky_residual_inf',resinf,flush=True)
print('C_certified_lower',clower,flush=True)
print('BLOCK_CERTIFIED',BLOCK,flush=True)
sys.exit(0)
