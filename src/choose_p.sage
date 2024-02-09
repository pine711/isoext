# This is an implementation to find the "good" prime for new method to compute kernel polynomials

from sage.all import *
#Given the embedding degree k, return candidates l
def candidates(k):
	res = []
	lowest = 2*k+1
	highest = floor(sqrt(k**3))
	for ii in range(lowest,highest):
		tmp = next_prime(ii)
		if k.divides(tmp-1):
			if tmp not in res:
				res.append(tmp)
	return res
def find_smallest_integer(l, k):
    """
    Input l and k, output the smallest integers a and b such that a, b are in the same coset
    with gcd(a,b)=1 and ord(a)|(l-1)/k.
    """
    F = GF(l)
    bound = min(l-1,k)
    A = []
    a0 = F.primitive_element()
    e = floor((l-1)/k)
    lam = a0**e
    S = []
    for x in range(0,k):
        s = lam**x
        S.append(s)
    for j in range(2,bound):
        a = F(j)
        orn_a = a.multiplicative_order()
        if orn_a % e ==0:
            B = []
            for x in range(0,k):
                b = S[x]*a
                B.append(b)
            flag = True
            while flag:
                b0 = min(B)
                b1 = max(B)
                if b0 == 1:
                    B.remove(b0)
                    b0 = min(B)
                if l-b1 == 1:
                    B.remove(b1)
                    b1 = max(B)
                if gcd(Integer(a),Integer(b0))==1:
                    flag = False
                    if b0 < a:
                        A.append((Integer(a),Integer(b0)))
                elif gcd(Integer(a),Integer(l-b1))==1:
                    flag = False
                    if Integer(l-b1)<Integer(a):
                        A.append((Integer(a),-Integer(l-b1)))
                else:
                    B.remove(b0)
    if len(A)!=0:
        (a,b)=min(A)
    else:
        a = 0
        b = 0
    return (a,b)

k = 332
l = 1993
k = 346
l = 3461
k = 470
l = 5641
k = 504
l = 2017
k = 544
l = 6529
# We want to implement this k and l.
k = 940
l = 28201
L = 14
n = 241 # the bit length of r, such that p = 2^L*r+1

find = False
iter = 0
print("ell = ", l, ", k = ", k, ", find p.")
while find == False:
	iter += 1
	r = ZZ(randint(2^(n-1),2^n-1))
	bound = 2^n-1-2^(n-1)
	p = r*2**L+1 # This case will make p = 1 mod 4, set p = r*2^L-1 can find p = 3 mod 4.
	if not p.is_pseudoprime():
		continue
	if p % 3 == 1 and p % 4 == 1:
		continue
	x = Mod(p, l).multiplicative_order() # p mod l, ord(p) in Z l
	if x%2 == 0:
		x //= 2
	if x == k:
		find = True
	if iter >= bound:
		print("Not find!")
		break
if find:
	print("Find!")
	print("Prime p = ",p)

# Find candidates k and l, such that max(|a|,|b|)=3.
init_k = 900
k = init_k
while k <= 1000:
    ls = candidates(k)
    for ii in range(len(ls)):
        a,b = find_smallest_integer(ls[ii],k)
        if max(abs(a),abs(b))==3:
            print("ell = ", ls[ii], ", k = ", k)
    k = k + 2