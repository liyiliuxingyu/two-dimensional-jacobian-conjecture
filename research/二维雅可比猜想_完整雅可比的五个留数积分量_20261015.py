#!/usr/bin/env python3
r"""
二维雅可比猜想续研：完整Jacobian的五个留数积分量
2026-10-15。研究工作记录，不是新排除论文、不是v8.57，不生成PDF/ZIP。

STATUS / 初学者入口
上轮是首项根的必要预筛。本轮尝试把所有剩余Jacobian方程一起处理，
而不是继续增加某一个系数或根例。所得是一个标准形式留数方法的等价改写，
不是新的排除结论、不主张方法或结论的优先性。真正cube/whole9/6/JC未解。

范围
K=C(t)，Q属于K[v]，deg_v Q=6，q6=H^6，H属于K*。
取s=Q^(1/6)~Hv，形式反演v=V(t,s)~s/H。
P=sum(j=0..9) c_j F_j，F_j=[Q^(j/6)]_+，c_j属于C，c9!=0。
保留全部十项，尤其c8；不套noncube的四项/反射/I方程。
由已恢复的完整Magnus接口，R=J(P,Q)的完整v次数<=4。
以下等价命题在K[v]成立；GLOBAL C[t,v]和原node是另加的必要条件。
有理/形式换元不是source多项式自同构，原52项的实际几何前提未取得。

精确结论（等价改写，不是“新排除”）
定义Res_v为形式Laurent级数的[v^-1]系数；不是几何无穷远留数的负号约定。
对m=1,...,5，定义
  M_m=(1/6)Res_v P Q^(m/6-1) Q_v
      =-(1/m)Res_v P_v Q^(m/6)。
它们属于K，是Phi(t,s)=P(t,V(t,s))的s^(-m)系数。
在上述完整候选形式下，对任意固定c属于C，有严格等价
  J(P,Q)=c  iff  M_1'=M_2'=M_3'=M_4'=0 且 M_5'=c/(6H)。
若要Keller搭配仍须c!=0；这一步不能省略。
在K=C(t)中，前四式等价于M_1,...,M_4都是复常数。
H=gamma(t-a)^k、k>=2时最后一式等价于
  M_5=d_5+c/[6 gamma(1-k)]*(t-a)^(1-k)，d_5属于C。
因此最后一个积分量有精确k-1阶极点；不能错误地要求五个都无极点/都是常数。
H常数时M_5=d_5+c*t/(6H)；k=1且c!=0时导数留数矛盾，不能有有理M_5。
这只重获既有1/H有理原函数门槛，不当作新单根分类。

PROOF / 完整证明
形式换元的留数不变性（首项非零线性）给
 [s^-m]Phi=Res_s Phi*s^(m-1) ds
 =Res_v P*Q^((m-1)/6)*(Q_v/(6Q^(5/6))) dv。
对v分部积分Res_v partial_v(P Q^(m/6))=0，得第二个表达。
固定s求t导数：V_t=-Q_t/Q_v，V_s=6s^5/Q_v，故
 Phi_t=R(t,V)/Q_v(t,V)=R(t,V)*V_s/(6s^5)。
再换元得完整恒等（对任意多项式P都成立，不先假R常数）
  M_m'=(1/6)Res_v R Q^(m/6-1)。
写R=sum(i=0..4) r_i v^i、Q=H^6 v^6 (1+E)，
 E=sum(i=0..5)(q_i/H^6)v^(i-6)。令A_alpha(n)=[z^n](1+E(z))^alpha，
 E(z)=sum(i=0..5)(q_i/H^6)z^(6-i)，A_alpha(n<0)=0。
则五个导数与(r0,...,r4)之间的矩阵为
  T_(m,i)=H^(m-6)/6*A_(m/6-1)(i+m-5)。
该矩阵反三角：m=1只含r4；m=2含r3/r4；依次类推。
反对角为H^-5/6,H^-4/6,H^-3/6,H^-2/6,H^-1/6，
行列式=H^-15/6^5（反序5阶的符号为正），所以没有漏掉任何r_i。
前四导数为0依次强迫r4=r3=r2=r1=0；第五式再强迫r0=c。
反向若R=c，前四留数次数小于-1而为0，第五恰为c/(6H)。等价得证。
K的导数常数域为C，故可积分为上述常数/显式单根表达。

重要适用性反例
deg_v R<=4不可删，单看五个M不能验证任意P的Jacobian。
例如Q=v^6+t，P=v^6-v，完整J=1-6v^5，并非常数；
Phi=s^6-t-(s^6-t)^(1/6)，所以M1..4=0、M5=t/6。
五个积分条件却与c=1相符。此P不满足十项常数组合接口（常数项含-t），
故不是本等价命题的反例，而是“不能跳过完整候选形式/degree界”的保护例。
其deg_v P=6，不冒充9/6反例、Keller模型或TEST。

全系数计算，不是增加喷射阶数
M_(j,m)=-(1/m)sum(i=1..j) i*[v^i]F_j * H^m * A_(m/6)(i+m)。
j<=9,m<=5时只需i+m<=14，这是形式留数的完整系数提取界，
不是以14阶近似代替未知更高项。j=0和j=6的全部M都恒0（P为1或Q）。
T的系数只用n<=4。对线性组合P，各M线性组合即可，没有未知Q统一次数界。

程序审计范围
脚本读取旧工具内固定SHA universal.json（该JSON已有不互导的二项式/六次根递推审计），
不是声称本脚本又独立证明一次Magnus一般定理。
这里另用(1+E)A'=alpha E'A的递推计算分数幂系数，
通用14自由变量H,q0..5及其独立一阶速度，标准库Fraction、无取样/浮点。
对全部50个M_(j,m)精确核验：其t导数=旧完整J_j经T变换；
核T的完整行列式、m=5的索引/符号、F0/F6零量及上述保护例。
复跑输出sha只证明字节稳定，不证明新排除或全未知Q有限。

NEXT / 未成功的尝试
把首项根高重数直接转为矛盾的尝试尚未完成。五个积分量只是把完整Jacobian
等价积分一次；前四常数和第五允许极点并不自动与全局多项式性冲突。
尤其不能把M5的必需极点误认为P或Q本身有极点。
下一步必须把这些精确有理恒等与全部P/Q多项式性、原node/完整support联合，
或取得actualF/H_B(F)/sourceCLCH/actual(e,f,k)/共同origintransport/源∞接触。
不要仅多算M、换记号、扩大保护例或重复等价接口升论文版本。
本轮没有新的通过Q、整族排除、原56/32/52删除；whole9/6/JC未解。

接续必须恢复的基线
最新首项根记录libfile_e8b4f4f9dc208191b5bff60c7706201f，
file_000000006840820ab7960bbd0cf19464，SHA b5b4ba62858d7a89ca80e54a0cca49e0a97a3ba1868dbd26eba0a28b408c4d83。
完整工具必须用修复v1：libfile_23d08a9a75f08191940c741e313fe986，
file_0000000034c481f5af73e5eaf5bd96cd，SHA e786b72f1f522dbe2f6a77165b514c777156bc190180a70d0d4fcf88c0c4cdcd。
本轮均已Library重新恢复并读完整记录、STATUS/NEXT/PROOF和嵌套v856 main。
v0坏ZIP已公开勘误，不再作为基线。最新排除论文仍v856，所有历史适用域保持。
2026-10-15重开Lee-Li arXiv2408.01279v1(2024-08-02)：
 https://arxiv.org/pdf/2408.01279 ，这是既有Magnus方法的背景，
本记录的留数换元/分部积分证明如上，不据其10.4单项式前提覆盖一般H。
本轮检索未取得可直接关闭当前GLOBAL9/6的已核primary定理；
搜索出现的网上突破或非学术站点声明不作为证明或状态输入。

运行：python3 本脚本.py /path/to/universal.json [--full]
默认输出全50项验证的摘要；--full还输出全部M和T的有理Laurent系数。
除输入文件外无依赖，不导入旧solver、不修改基线。
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import permutations
import json,sys,hashlib
N=14
ZERO=(0,)*N
def add(*ps):
    r={}
    for p in ps:
        for a,c in p.items():r[a]=r.get(a,F(0))+c
    return {a:c for a,c in r.items() if c}
def sc(p,c):return {a:c*d for a,d in p.items() if c*d}
def mul(p,q):
    r={}
    for a,c in p.items():
        for b,d in q.items():
            k=tuple(x+y for x,y in zip(a,b));r[k]=r.get(k,F(0))+c*d
    return {a:c for a,c in r.items() if c}
def shift(p,h):
    return {(a[0]+h,)+a[1:]:c for a,c in p.items()}
def var(i):
    a=[0]*N;a[i]=1;return {tuple(a):F(1)}
def dt(p):
    r={}
    for a,c in p.items():
        assert not any(a[7:]),'input must have no velocity variables'
        for i in range(7):
            if a[i]:
                b=list(a);b[i]-=1;b[i+7]+=1;k=tuple(b);r[k]=r.get(k,F(0))+c*a[i]
    return {a:c for a,c in r.items() if c}
def decode(terms):
    r={}
    for ex,cs in terms:r.setdefault(ex[14],{})[tuple(ex[:14])]=F(cs)
    return r
def series(alpha,L):
    # z^i coefficient of E is q_(6-i)/H^6, i=1,...,6.
    e={i:shift(var(7-i),-6) for i in range(1,7)}
    a=[{ZERO:F(1)}]
    for n in range(1,L+1):
        a.append(sc(add(*[sc(mul(e[i],a[n-i]),(alpha+1)*i-n) for i in range(1,min(6,n)+1)]),F(1,n)))
    return a
def encode(p):return [[list(a),str(c)] for a,c in sorted(p.items())]
def verify(path,full=False):
    data=Path(path).read_bytes()
    assert hashlib.sha256(data).hexdigest()=='b5977dda6b6c8f1c009abec8b5a189bd4d7d69daefe08bca71c4745943c9910a'
    u=json.loads(data);basis=[decode(x) for x in u['basis']];js=[decode(x) for x in u['jacobians']]
    assert all(all(i<=4 for i in p) for p in js)
    T=[];moments={}
    for m in range(1,6):
        a=series(F(m,6),14);b=series(F(m,6)-1,4)
        row=[sc(shift(b[i+m-5],m-6),F(1,6)) if i+m>=5 else {} for i in range(5)]
        T.append(row)
        for j in range(10):
            moment=sc(shift(add(*[sc(mul(p,a[i+m]),i) for i,p in basis[j].items() if i]),m),F(-1,m))
            rhs=add(*[mul(row[i],p) for i,p in js[j].items()])
            assert dt(moment)==rhs,(j,m)
            if j in (0,6):assert not moment
            moments[f'F{j}_M{m}']=moment
    det={}
    for perm in permutations(range(5)):
        z={ZERO:F(1)}
        for i,j in enumerate(perm):z=mul(z,T[i][j])
        sign=(-1)**sum(perm[i]>perm[j] for i in range(5) for j in range(i+1,5))
        det=add(det,sc(z,sign))
    assert det=={(-15,)+(0,)*13:F(1,6**5)}
    for m,row in enumerate(T,1):
        assert row[0]==({(-1,)+(0,)*13:F(1,6)} if m==5 else {})
    # General R=c+v^5 lies outside degree<=4 and its v^5 term is invisible
    # to all five moments for Q=v^6+t. The explicit P=v^6-v has c=1.
    # Verify directly with univariate integer exponent dictionaries in t,v.
    P={(0,6):F(1),(0,1):F(-1)};Q={(0,6):F(1),(1,0):F(1)}
    def d(p,k):
        r={}
        for a,c in p.items():
            if a[k]:b=list(a);b[k]-=1;r[tuple(b)]=c*a[k]
        return r
    def product(p,q):
        r={}
        for a,c in p.items():
            for b,e in q.items():k=(a[0]+b[0],a[1]+b[1]);r[k]=r.get(k,F(0))+c*e
        return {a:c for a,c in r.items() if c}
    raw=product(d(P,0),d(Q,1))
    for a,c in product(d(P,1),d(Q,0)).items():raw[a]=raw.get(a,F(0))-c
    assert {a:c for a,c in raw.items() if c}=={(0,0):F(1),(0,5):F(-6)}
    for m in range(1,6):
        # Directly use -Res(P_v Q^(m/6))/m, not the asserted Phi formula.
        found={}
        for vi,coeff in [(5,F(6)),(0,F(-1))]:
            if (vi+m+1)%6:continue
            n=(vi+m+1)//6;b=F(1)
            for k in range(n):b*=F(F(m,6)-k,k+1)
            found[n]=found.get(n,F(0))-coeff*b/m
        found={n:c for n,c in found.items() if c}
        assert found==({1:F(1,6)} if m==5 else {})
    full_obj=dict(moments={k:encode(p) for k,p in sorted(moments.items())},matrix=[[encode(p) for p in row] for row in T],determinant=encode(det))
    canonical=json.dumps(full_obj,sort_keys=True,separators=(',',':')).encode()
    out=dict(status='PASS',full_moment_identities=50,non_sample_symbolic=True,
             determinant='H^-15/6^5',coefficients_sha256=hashlib.sha256(canonical).hexdigest(),
             degree_bound_guard='P=v^6-v,Q=v^6+t: moments pass c=1 but J=1-6v^5',
             new_exclusion_claimed=False)
    if full:out.update(full_obj)
    return out
if __name__=='__main__':
    if len(sys.argv) not in (2,3) or (len(sys.argv)==3 and sys.argv[2]!='--full'):
        raise SystemExit('usage: python3 script.py universal.json [--full]')
    print(json.dumps(verify(sys.argv[1],len(sys.argv)==3),ensure_ascii=False,sort_keys=True,indent=2))
