from deuring.broker import starting_curve
from deuring.randomideal import random_ideal
from deuring.correspondence import constructive_deuring, constructive_deuring_new
from sage.all import *

x = var("x")
#p3923
# p = 23759399264157352358673788613307970528646815114090876784643387662192449945599
# F2 = GF(p**2, name = "i", modulus = x**2+1)
#p2
p = 37670568336551536389503919665937491111216122470333837677213877442445311999999
F2 = GF(p**2, name = "i", modulus = x**2+1)
#p1
# p = 11956566944641502957704189594909498993478297403838643406058180334130656750161
# F2 = GF(p**2, name = "i", modulus = x**2+17)

t = cputime()
E0, iota, O0 = starting_curve(F2)
print("The runtime of constructing E0:",cputime(t),"s")
# print(E0)
# print(O0)
t1 = 0
t2 = 0
num = 1
for ii in range(num):
    I = random_ideal(O0)
    t = cputime()
    E2, phi2, _ = constructive_deuring_new(I, E0, iota)
    t1+=cputime(t)
    t = cputime()
    E1, phi, _ = constructive_deuring(I, E0, iota)
    t2+=cputime(t)
print("==================================================")
print("The original algorithem for computing constructive deuring correspondence takes",t2/num,"s")
print("The new algorithem for computing constructive deuring correspondence takes",t1/num,"s")