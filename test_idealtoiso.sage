from deuring.broker import starting_curve
from sage.all import *
import collections
from deuring.idealtoisogeny import Deuring_Context

x = var("x")
# p = 1 mod 4, p = 2 mod 3, k = 252, l = 1009
p = 49139181582706039787278142920262584477092436049503079043059897492654038958081
l = ZZ(1009)
k = 252
# p = 1 mod 4, p = 2 mod 3, k = 346, l = 3461
# p = 39080947418725606808429296766910149105785171601742635268655886999006895616001
# l = ZZ(3461)
# k = 346
# p = 1 mod 4, p = 2 mod 3, k = 470, l = 5641
# p = 51980146777764317376899789838214467900034016312224807258633370210069137088513
# l = ZZ(5641)
# k = 470
# p = 1 mod 4, p = 2 mod 3, k = 504, l = 2017
# p = 49514069746810538044618160680453043783501727017531641497308339189045991350273
# l = ZZ(2017)
# k = 504
# p = 1 mod 4, p = 2 mod 3, k = 940, l = 28201
# p = 56893651996855845314990170884559497875088198021359275302882959317982204772353
# l = ZZ(28201)
# k = 940
# p = 1 mod 4, p = 2 mod 3, k = 948, l = 3793
# p = 29270790387301979192886935629516082056361570481087354970797119294565959950337
# l = ZZ(3793)
# k = 948
F2 = GF(p**2, name = "i", modulus = x**2+3)
t = cputime()
E0, iota, O0 = starting_curve(F2)
print("The runtime of constructing E0:",cputime(t),"s")
print(E0)
print(O0)
################################
# Random a left ideal I of norm l.
while True:
    alpha = O0.random_element()
    if gcd(alpha.reduced_norm(),l)==l and alpha!=0:
        break
#alpha = -20*O0.gens()[0]-O0.gens()[2]+13*O0.gens()[3] #k=940
I = O0*alpha + O0*l
facToExt = dict()
facToExt[l] = k
print(f'norm(I) = {factor(I.norm())}')
grouped = collections.defaultdict(list)
for le,k in sorted(facToExt.items()):
    grouped[k].append(le)
print('Torsion by extension degree:')

for k,les in sorted(grouped.items()):
    print(f'  {k:4}:  ', end='')
    for j,le in enumerate(les):
        print((', ' if j else '') + str(factor(le)), end='')
    print()
Deuring_ctx = Deuring_Context(O0, E0, iota, facToExt)
# pari.allocatemem(23554432000000) # For big k, we may need to allocate more memory
phi_I = Deuring_ctx.IdealToIsogeny(I)