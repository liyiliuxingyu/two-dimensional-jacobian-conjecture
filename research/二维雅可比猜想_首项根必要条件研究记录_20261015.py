#!/usr/bin/env python3
r"""
二维雅可比猜想：真正立方首项根处的必要条件研究记录
日期：2026-10-15。不是v8.57论文，不新增整族排除，不生成PDF。

初学者入口 / STATUS
上一工具检查一个完整Q。本记录先问：不猜P、也不给未知Q设总次数界，
有哪些系数形状绝不可能通过？得到的是必要条件，不是存在判据。
以下均为既有Newton/Magnus方法的本层推论；不主张发现或优先性。
通过本脚本只表示未被这些条件拒绝，绝非Keller候选认证。

完整范围
P,Q属于C[t,v]，deg_v P=9，deg_v Q=6，q6=H^6，0!=H属于C[t]，
J(P,Q)=P_t Q_v-P_v Q_t属于C*。不需要node直线或t次数界。
保留全部十项P=sum c_j F_j，F_j=[Q^(j/6)]_+，c_j属于C，c9!=0。
特别保留c8；不能用非cube的四常数、反射或I方程代替。
CLI仅实现有理系数H/q0..q5的必要条件预筛，精确mate_degree必须9。

本轮工作结论（不是新的全层排除）
(A) 每个a属于C，(deg_v P(a,v),deg_v Q(a,v))属于{(3,2),(6,4),(9,6)}。
    H(a)=0时只剩前两项；这是已有Newton多边形相似定理的应用。
(B) H^3整除q5。置B=q5/H^3，则B是多项式。
(C) H整除E=q4-B^2/4。置C=E/H，则C是多项式。
(D) 必有同一个常数tau=c8/c9满足
        243 C^2+40 tau B^3 == 0 (mod H)。
    不是允许tau随t改变，也不能在各根独立选择tau。
(E) 令R=H的monic平方自由部，S=gcd(R,q4)。则
        S整除q3，gcd(S,q2)=1。
    (C)已保证q4在H的根上等于B^2/4，所以也可用S=gcd(R,B)。
    R/S上的纤维Q次数4，S上的纤维Q次数2。H常数时均为真空条件。

完整证明 / PROOF
一、(A)/(E)。用NP(f)=conv(supp(f) union {(0,0)})，不是不加原点的凸包。
已知全局常Jac、两个总次数>1时，NP(P)和NP(Q)以原点为中心相似。
最大v坐标9和6决定NP(P)=(3/2)NP(Q)，也决定总次数比3/2。
对任意a只作源仿射变化x=t-a，仍满足同一已知定理和v高度9/6。
由于x指数非负，NP与v轴的截距就是各限制P(a,v)、Q(a,v)的次数
（若限制为零或常数，截距记0）。两截距之比仍为3/2。
它们不能同时为0，否则P_v和Q_v在整条t=a线上为0，J也为0。
故截距只能是(3,2)、(6,4)、(9,6)。H(a)=0使Q次数<6，排除最后一项。
于是q5(a)=0；q4(a)=0时又须q3(a)=0、q2(a)!=0。这正是(E)。
源平移是全局多项式仿射变化，不是把有理中心当作多项式坐标。
此处“纤维”指原源直线t=a上的限制，不是目标Q的水平纤维。

二、(B)的独立允许项证明。非负主部的完整v6系数为
 [v6]F9=3H^3 q3/2+3 q5 q4/(4H^3)-q5^3/(16H^9),
 [v6]F8=4H^2 q4/3+2 q5^2/(9H^4),
 [v6]F7=7H q5/6，[v6]F6=H^6；j<6不贡献。
设有限根a处ord(H)=k>=1，q5非零且m=ord(q5)<3k。
F9的立方项阶3m-9k<0。其他可能含极点项的阶分别至少
 m-3k、2m-4k，其余项非负；因
 3m-9k < m-3k、3m-9k < 2m-4k，且c9!=0，最低阶唯一，
不可能抵消为多项式。故每个根m>=3k，得到H^3|q5。
q5恒0直接满足，没有将零支取对数或偷删。

三、(C)/(D)。代入q5=H^3 B、q4=B^2/4+E，完整v5系数为
 [v5]F9=3H^3 q2/2+3B q3/4+3E^2/(8H^3),
 [v5]F8=4H^2 q3/3+4BE/(9H)+5B^3/(81H),
 [v5]F7=H(7E/6+7B^2/18)，[v5]F6=H^3 B，[v5]F5=H^5。
j<5不贡献。若a处e=ord(E)<k=ord(H)，则c9项的阶2e-3k<-k，
c8项阶至少-k，其余为多项式。再次唯一最低极点，矛盾，故H|E。
E恒0亦直接保留。置E=H C后，P的v5系数唯一可能的分母部分为
       (243 c9 C^2+40 c8 B^3)/(648H)。
其余全部为多项式，所以在(B)/(C)之后，v5系数多项式性等价于(D)。
这一“等价”只针对该单个系数，绝非整个P多项式性或常Jac的充要条件。
仍需原工具的所有余式、node限制和完整Jac条件；更低允许项从未删除。

四、计算边界。有限个有理输入只有这里的常数tau被判定，没有枚举未知Q。
若243 C^2 mod H和40 B^3 mod H均为0，tau任意；
否则对后者任一非零系数求tau并核全部系数即可，余式不成比例则拒绝。
有理输入时若有复tau，此处也有有理tau（或tau任意）；不是只排有理常数。
主部系数恒等式用标准库Fraction的通用Laurent环逐项核验，不靠有限取样。
下面自测的通过输入只是必要条件的反误判保护，不是常Jac对或实际Q实现。

NEXT / 未解决
尚未证明(B)-(E)能迫使所有更低主部闭合，更未关闭真正cube九/六。
不要把这些必要条件、通过样例、tau存在或部分系数整除当作充分条件。
后续应联结其余完整多项式系数/Jac或真实源数据；不要仅扩大本自测表发论文。
完整TEST仍需v850支撑/k>=11；本记录不依赖该下界，不认证任何原CL可化TEST。
actualF/H_B(F)/sourceCLCH/actual(e,f,k)/共同origintransport/源∞接触仍缺，
原56/32/52无新增删除，whole9/6和JC未解。其余严格历史边界保持。

接续来源（精确基线另存，避免重复封装大包）
给定Q工具最终版本1：libfile_23d08a9a75f08191940c741e313fe986，
file_0000000034c481f5af73e5eaf5bd96cd，
SHA e786b72f1f522dbe2f6a77165b514c777156bc190180a70d0d4fcf88c0c4cdcd。
本轮已重新恢复完整CORRECTION/STATUS/NEXT/PROOF及嵌套v856 STATUS/NEXT/main；
只用修复版1，不用已损坏v0。最新排除论文仍v856，不重包旧结论。
2026-10-15重新打开primary：
Gwozdziewicz, Injectivity on one line (1993), Theorem2.1，全文两页，
 https://arxiv.org/pdf/alg-geom/9305008 。该文引Abhyankar/Lang原定理；
本记录不声称重证Lang或读过其原稿。
Lee-Li, Magnus' formula revisited IV, arXiv2408.01279v1 (2024-08-02)，
 https://arxiv.org/pdf/2408.01279 。允许项直接由已保存形式反演证明；
不借其10.4单项式条件覆盖任意H，不声称新Magnus一般定理。

运行：python3 本文件.py         # 通用恒等式+有限实现自测，输出JSON
      python3 本文件.py input.json  # 与旧工具同H/q数组；只做上述必要预筛
      python3 本文件.py --audit-baseline universal.json
        # 与旧工具内固定SHA、曾由另一套六次根递推得到的通用式比较
不导入旧求解器，不用浮点、有限域或根的数值近似。
"""
from fractions import Fraction as F
import json, sys, hashlib
from pathlib import Path

# Sparse Laurent ring in H,q0,...,q5,B,E, then coefficient in v.
DIM=9
def add(*ps):
    r={}
    for p in ps:
        for k,c in p.items():r[k]=r.get(k,F(0))+c
    return {k:c for k,c in r.items() if c}
def scale(p,c):return {k:c*x for k,x in p.items() if c*x}
def mul(p,q):
    r={}
    for a,c in p.items():
        for b,d in q.items():
            k=tuple(x+y for x,y in zip(a,b));r[k]=r.get(k,F(0))+c*d
    return {k:c for k,c in r.items() if c}
def power(p,n):
    r={(0,)*DIM:F(1)}
    for _ in range(n):r=mul(r,p)
    return r
def var(i,n=1):
    a=[0]*DIM;a[i]=n;return {tuple(a):F(1)}
def h(n):return var(0,n)
def q(i):return var(i+1)
def binomial(x,n):
    r=F(1)
    for i in range(n):r*=F(x-i,i+1)
    return r
def generic(j,vi):
    drop=j-vi
    if drop<0:return {}
    terms={6-i:q(i) for i in range(6)};ps={0:{(0,)*DIM:F(1)}};result={}
    for n in range(drop+1):
        result=add(result,scale(mul(h(j-6*n),ps.get(drop,{})),binomial(F(j,6),n)))
        nxt={}
        for a,p in ps.items():
            for b,r in terms.items():
                if a+b<=drop:nxt[a+b]=add(nxt.get(a+b,{}),mul(p,r))
        ps=nxt
    return result
def subst(p):
    values=[var(i) for i in range(DIM)]
    values[6]=mul(h(3),var(7))
    values[5]=add(scale(power(var(7),2),F(1,4)),var(8))
    r={}
    for ex,c in p.items():
        z={(0,)*DIM:c}
        for i,n in enumerate(ex):
            z=mul(z,h(n) if i==0 else power(values[i],n))
        r=add(r,z)
    return r
def identities():
    B,E=var(7),var(8)
    literal6={
      9:add(scale(mul(h(3),q(3)),F(3,2)),scale(mul(h(-3),mul(q(5),q(4))),F(3,4)),scale(mul(h(-9),power(q(5),3)),F(-1,16))),
      8:add(scale(mul(h(2),q(4)),F(4,3)),scale(mul(h(-4),power(q(5),2)),F(2,9))),
      7:scale(mul(h(1),q(5)),F(7,6)),6:h(6)}
    literal5={
      9:add(scale(mul(h(3),q(2)),F(3,2)),scale(mul(B,q(3)),F(3,4)),scale(mul(h(-3),power(E,2)),F(3,8))),
      8:add(scale(mul(h(2),q(3)),F(4,3)),scale(mul(h(-1),mul(B,E)),F(4,9)),scale(mul(h(-1),power(B,3)),F(5,81))),
      7:mul(h(1),add(scale(E,F(7,6)),scale(power(B,2),F(7,18)))),6:mul(h(3),B),5:h(5)}
    rows={}
    for vi,table in [(6,literal6),(5,literal5)]:
        for j in range(10):
            got=generic(j,vi)
            if vi==5:got=subst(got)
            assert got==table.get(j,{})
            rows[f'F{j}_v{vi}']=[[list(k),str(c)] for k,c in sorted(got.items())]
    return rows

# Ordinary Q[t] arithmetic; no root factorization or rational-function black box.
def trim(p):
    while p and not p[-1]:p.pop()
    return p
def plus(p,q):
    r=[F(0)]*max(len(p),len(q))
    for i,c in enumerate(p):r[i]+=c
    for i,c in enumerate(q):r[i]+=c
    return trim(r)
def times(p,q):
    r=[F(0)]*max(0,len(p)+len(q)-1)
    for i,c in enumerate(p):
        for j,d in enumerate(q):r[i+j]+=c*d
    return trim(r)
def ppow(p,n):
    r=[F(1)]
    for _ in range(n):r=times(r,p)
    return r
def divrem(p,q):
    if not q:raise ValueError('division by zero')
    p=p[:];r=[F(0)]*max(0,len(p)-len(q)+1)
    while p and len(p)>=len(q):
        k=len(p)-len(q);c=p[-1]/q[-1];r[k]=c
        for j,b in enumerate(q):p[k+j]-=c*b
        trim(p)
    return trim(r),p
def monic(p):return [c/p[-1] for c in p] if p else []
def gcd(p,q):
    while q:p,q=q,divrem(p,q)[1]
    return monic(p)
def derivative(p):return trim([i*p[i] for i in range(1,len(p))])
def encoded(p):return [str(c) for c in p]
def parse(p):
    if not isinstance(p,list) or any(type(c) not in (str,int) for c in p):raise ValueError('exact integer/rational arrays required')
    return trim([F(c) for c in p])
def preflight(inp):
    N=inp.get('mate_degree',9)
    if type(N) is not int or N!=9:raise ValueError('this necessary preflight is only for exact mate degree 9')
    H=parse(inp['H']);qs=inp['q']
    if not H or not isinstance(qs,list) or len(qs)!=6:raise ValueError('nonzero H and exactly q0..q5 required')
    qs=[parse(p) for p in qs]
    out={'scope':'necessary conditions only; not closure or existence','node_constraint_used':False}
    B,r=divrem(qs[5],ppow(H,3))
    if r:return dict(out,necessary_pass=False,failed='H^3 divides q5',remainder=encoded(r))
    E=plus(qs[4],[-c/4 for c in ppow(B,2)])
    C,r=divrem(E,H)
    if r:return dict(out,necessary_pass=False,failed='H divides q4-B^2/4',remainder=encoded(r))
    r=divrem([243*c for c in ppow(C,2)],H)[1]
    s=divrem([40*c for c in ppow(B,3)],H)[1]
    if not s:
        if r:return dict(out,necessary_pass=False,failed='single constant tau congruence',R=encoded(r),S=encoded(s))
        tau='any'
    else:
        i=next(i for i,c in enumerate(s) if c);tau=-(r[i] if i<len(r) else F(0))/s[i]
        if plus(r,[tau*c for c in s]):return dict(out,necessary_pass=False,failed='single constant tau congruence',R=encoded(r),S=encoded(s))
        tau=str(tau)
    R=monic(divrem(H,gcd(H,derivative(H)))[0]);S=gcd(R,qs[4])
    rem=divrem(qs[3],S)[1]
    if rem:return dict(out,necessary_pass=False,failed='fiber degree 3 forbidden',fiber_factor=encoded(S),remainder=encoded(rem))
    common=gcd(S,qs[2])
    if len(common)>1:return dict(out,necessary_pass=False,failed='fiber degree at most 1 forbidden',fiber_factor=encoded(common))
    assert gcd(R,B)==S
    return dict(out,necessary_pass=True,B=encoded(B),C=encoded(C),tau=tau,
                fiber_degree2_factor=encoded(S),fiber_degree4_factor=encoded(divrem(R,S)[0]))
def selftest():
    def case(H,q2=None,q3=None,q4=None,q5=None):return dict(H=H,q=[[],[],q2 or [],q3 or [],q4 or [],q5 or []],mate_degree=9)
    tests=[
      (case([0,1],q2=[1]),True),
      (case([0,1],q4=[1],q5=[0,0,0,2]),True),
      (case([0,1],q5=[1]),False),
      (case([0,0,1],q4=[0,1]),False),
      (case([0,0,1],q4=[F(1,4),0,1,1],q5=[0,0,0,0,0,0,1]),False),
      (case([0,1],q2=[1],q3=[1]),False),
      (case([0,1]),False)]
    for inp,expected in tests:
        inp['q']=[[str(c) for c in p] for p in inp['q']]
        assert preflight(inp)['necessary_pass']==expected
    # Floating point inputs must not enter an exact certificate.
    try:preflight(case([0.0,1],q2=[1]))
    except ValueError:pass
    else:raise AssertionError('accepted float')
    return {'passed':len(tests),'float_rejected':True,'positive_inputs_are_not_Keller_pairs':True}
def audit_baseline(path):
    data=Path(path).read_bytes()
    assert hashlib.sha256(data).hexdigest()=='b5977dda6b6c8f1c009abec8b5a189bd4d7d69daefe08bca71c4745943c9910a'
    old=json.loads(data)
    for j in range(10):
        for vi in (5,6):
            p={}
            for ex,cs in old['basis'][j]:
                if ex[14]!=vi:continue
                assert all(e==0 for e in ex[7:14])
                p[tuple(ex[:7])+ (0,0)]=F(cs)
            assert p==generic(j,vi)
    return {'twenty_generic_coefficients_match_pinned_baseline':True,
            'baseline_is_not_a_new_independent_derivation':True}
if __name__=='__main__':
    rows=identities()
    if len(sys.argv)==1:
        out=dict(generic_identities=rows,implementation_tests=selftest(),status='PASS',new_exclusion_claimed=False)
    elif len(sys.argv)==2:
        out=preflight(json.loads(Path(sys.argv[1]).read_text()))
    elif len(sys.argv)==3 and sys.argv[1]=='--audit-baseline':
        out=audit_baseline(sys.argv[2]);out.update(status='PASS',implementation_tests=selftest())
    else:raise SystemExit('usage: python3 script.py [input.json]')
    print(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2))
