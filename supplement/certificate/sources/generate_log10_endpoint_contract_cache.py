"""Generate fixed log10/log11 endpoint contract caches; scratch gate only.

The historical positional interface is unchanged.  Pass ``--endpoint 11``
to select ``a=log(11)/2`` and its independently certified prime bound.
"""
import sys

if '--endpoint-rational' in sys.argv:
    ri=sys.argv.index('--endpoint-rational');ENDPOINT_NUM=int(sys.argv[ri+1]);ENDPOINT_DEN=int(sys.argv[ri+2])
    if ENDPOINT_NUM <= 0 or ENDPOINT_DEN <= 0:
        raise ValueError('rational endpoint must be positive')
    ENDPOINT=None;CACHE_PREFIX=f'log{ENDPOINT_NUM}over{ENDPOINT_DEN}'
    PRIME_BOUND='2.931158';A_EXPR=f'(arb({ENDPOINT_NUM})/arb({ENDPOINT_DEN})).log()/2'
else:
    ENDPOINT = (int(sys.argv[sys.argv.index('--endpoint') + 1])
                if '--endpoint' in sys.argv else 10)
    if ENDPOINT not in {10, 11}: raise ValueError('supported endpoints are 10 and 11')
    CACHE_PREFIX = f'log{ENDPOINT}';PRIME_BOUND = {10: '2.978071', 11: '2.931158'}[ENDPOINT]
    A_EXPR=f'arb({ENDPOINT}).log()/2'
if '--prime-bound' in sys.argv:
    pi=sys.argv.index('--prime-bound')
    PRIME_BOUND=sys.argv[pi+1]
TAIL_GAP='1e-12'
if '--tail-gap' in sys.argv:
    ti=sys.argv.index('--tail-gap')
    if ti+1 >= len(sys.argv):
        raise ValueError('--tail-gap requires a positive decimal')
    TAIL_GAP=sys.argv[ti+1]
    if not TAIL_GAP or TAIL_GAP.startswith('-'):
        raise ValueError('--tail-gap must be positive')
print('tail_gap', TAIL_GAP, flush=True)
GRAPH_HOOK = f"""if '--export-verified-graph' in sys.argv and R==280 and N==500:
        tag='even' if PAR=='+' else 'odd'
        np.savez(f'scratch/{CACHE_PREFIX}_{{tag}}_R280_N500_verified_graph_arb_{{ctx.prec}}.npz',Z0=np.array([[str(zcols[j][i]) for j in range(BLOCK)] for i in range(R)]),S=np.array([[str(small[i,j]) for j in range(BLOCK)] for i in range(BLOCK)]),X=np.array([[str(X[i,j]) for j in range(m)] for i in range(BLOCK)]),C=np.array([[str(Cmat[i,j]) for j in range(m)] for i in range(m)]),whist=np.array([[str(w[i]) for i in range(R)] for w in whist]),precision=np.array(ctx.prec),parity=np.array(PAR),R=np.array(R),N=np.array(N),a=np.array(str(a)),mu=np.array(str(mu)),cache_prefix=np.array(CACHE_PREFIX),prime_bound=np.array(str(prime_bound)))
        print('EXPORTED_VERIFIED_GRAPH',tag,ctx.prec,flush=True);raise SystemExit(0)
    if '--export-small-map-only' in sys.argv:"""
BODY_CODE = (
    f"CACHE_PREFIX={CACHE_PREFIX!r};"
    "body=open('computations/log8_block_certificate_body.py',encoding='utf-8').read();"
    "body=body.replace('assert PAR == \"+\" and R == 128','assert R == 280').replace(\"assert PAR == '+' and R == 128\",'assert R == 280')"
    ".replace('BLOCK = 4','BLOCK = (7 if PAR==\"+\" else 6)')"
    ".replace('CHOLESKY_BITS = 52','CHOLESKY_BITS = 96')"
    f".replace('C0 = arb(1) / 2000000',\"C0 = arb({TAIL_GAP!r})\")"
    ".replace('assert clower.lower()>0','print(\"precheck_resinf\",resinf,\"precheck_clower\",clower,flush=True);assert clower.lower()>0')"
    ".replace('if export_case and BLOCK==5 and N==400:','if N==500:')"
    ".replace(\"if PAR=='-' and R==200 and BLOCK==5:\",\"if False:\")"
    f".replace(\"if '--export-small-map-only' in sys.argv:\",{GRAPH_HOOK!r});"
    "exec(compile(body,__file__,'exec'))"
)
with open('computations/verify_first_mixed_band_endpoint_log5.py',encoding='utf-8') as h:src=h.read()
changes={
'pi=arb.pi(); a=arb(5).log()/2':f'pi=arb.pi(); a={A_EXPR}',
    'd2=arb(2).log();d3=arb(3).log();d4=2*d2':'d2=arb(2).log();d3=arb(3).log();d4=2*d2;d5=arb(5).log();d7=arb(7).log();d8=3*d2;d9=2*d3;d11=arb(11).log()',
    'alpha2=d2/arb(2).sqrt();alpha3=d3/arb(3).sqrt();alpha4=d2/2':'alpha2=d2/arb(2).sqrt();alpha3=d3/arb(3).sqrt();alpha4=d2/2;alpha5=d5/arb(5).sqrt();alpha7=d7/arb(7).sqrt();alpha8=d2/arb(8).sqrt();alpha9=d3/3;alpha11=d11/arb(11).sqrt()',
    "TERMS += [(-eps*alpha3,'ref',d3,d3-a,a),(-eps*alpha4,'ref',d4,d4-a,a)]":"TERMS += [(-eps*alpha3,'ref',d3,arb(0),d3),(-alpha3,'plus',d3,arb(0),a-d3),(-alpha3,'minus',d3,d3,a),(-eps*alpha4,'ref',d4,d4-a,a),(-eps*alpha5,'ref',d5,d5-a,a),(-eps*alpha7,'ref',d7,d7-a,a),(-eps*alpha8,'ref',d8,d8-a,a),(-eps*alpha9,'ref',d9,d9-a,a),(-eps*alpha11,'ref',d11,d11-a,a)]",
'prime_bound=2*alpha2+alpha3+alpha4':f'prime_bound=arb("{PRIME_BOUND}");print("prime_bound",prime_bound,flush=True)',
"L=[[arb(0) for _ in range(R)] for _ in range(R)];dd=[arb(0) for _ in range(R)];mn=None":BODY_CODE,
}
if '--export-full-safe-matrix' in __import__('sys').argv:
 key="L=[[arb(0) for _ in range(R)] for _ in range(R)];dd=[arb(0) for _ in range(R)];mn=None"
 changes[key]=f"numpy=__import__('numpy');numpy.savez(f'scratch/{CACHE_PREFIX}_{{\"even\" if PAR==\"+\" else \"odd\"}}_K{{R}}_safe_base_arb_{{ctx.prec}}.npz',S=numpy.array([[str(S[i,j]) for j in range(R)] for i in range(R)]),Rcov=numpy.array([[str(CC[i,j]) for j in range(R)] for i in range(R)]),parity=numpy.array(PAR),R=numpy.array(R),N=numpy.array(N),precision=numpy.array(ctx.prec),a=numpy.array(str(a)),mu=numpy.array(str(mu)),prime_bound=numpy.array(str(prime_bound)),cache_prefix=numpy.array(CACHE_PREFIX));print('EXPORTED_FULL_SAFE_MATRIX',R,ctx.prec,flush=True);raise SystemExit(0)"
if '--evaluate-dyadic-witness' in __import__('sys').argv:
 key="L=[[arb(0) for _ in range(R)] for _ in range(R)];dd=[arb(0) for _ in range(R)];mn=None"
 hook="exec(compile(open('scratch/export_dyadic_subspace_quadratics.py',encoding='utf-8').read(),__file__,'exec'))"
 changes[key]=hook
if '--export-r500-parseval' in __import__('sys').argv:
 key="L=[[arb(0) for _ in range(R)] for _ in range(R)];dd=[arb(0) for _ in range(R)];mn=None"
 changes[key]="exec(compile(open('scratch/export_log10_even_r500_parseval_data.py',encoding='utf-8').read(),__file__,'exec'))"
if '--export-yukawa-family' in __import__('sys').argv:
 argv=__import__('sys').argv;flag=argv.index('--export-yukawa-family')
 yukawa_lo=int(argv[flag+1]);yukawa_hi=int(argv[flag+2])
 assert 500 <= yukawa_lo < yukawa_hi
 key="L=[[arb(0) for _ in range(R)] for _ in range(R)];dd=[arb(0) for _ in range(R)];mn=None"
 changes[key]=(f"YUKAWA_EXPORT_LO={yukawa_lo};YUKAWA_EXPORT_HI={yukawa_hi};"
               "exec(compile(open('scratch/export_log10_even_yukawa_family_parseval_data.py',encoding='utf-8').read(),__file__,'exec'))")
for old,new in changes.items():
 if "body=open('computations/log8_block_certificate_body.py'" in new:
  new=new.replace(".replace('if export_case and BLOCK==5 and N==400:'", ".replace(\"f'_R{R}_N400_exact_small_map.npz'\",\"f'_R{R}_N{N}_exact_small_map.npz'\").replace('if export_case and BLOCK==5 and N==400:'")
 assert src.count(old)==1,(old,src.count(old));src=src.replace(old,new)
exec(compile(src,__file__,'exec'))
