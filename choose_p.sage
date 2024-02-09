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
    with gcd(a,b)=1 and (l-1)/k|ord(a).
    """
    assert k>1
    F = GF(l)
    bound = l # Find a and b in Z_\ell
    A = []
    # Find a primitive root of Z_\ell
    from sage.schemes.elliptic_curves.isogeny_small_degree import _least_semi_primitive
    a0 = F(_least_semi_primitive(l))
    e = (l-1)//k
    # Construct the set S = {1,\lambda,\lambda^2,...,\lambda^{k-1}}
    lam = a0**e
    S = []
    for x in range(0,k):
        s = lam**x
        S.append(s)
    # Find a and b
    R = []
    for j in range(2,bound):
        a = F(j)
        ord_a = a.multiplicative_order()
        if ord_a % e ==0:
            A = []
            B = []
            for x in S:
                b = x*a
                B.append(b)
            # Filter B by <+->
            B1 = []
            for x in B:
                tmp = (l-1)//2
                if x > tmp:
                    B1.append(Integer(l-x))
                else:
                    B1.append(Integer(x))
            B1.sort() # Sort the set of B1 with ascending order
            B1 = set(B1) # Filter repeated elements
            if 1 in B1:
                B1.remove(1)
            while len(B1)!=0:
                candidate = B1.pop()
                if gcd(Integer(a),candidate) == 1:
                    A.append(candidate)
                    break
        if len(A)!=0:
            R.append((a,A[0]))
    if len(R)!=0:
        a,b = min(R)
    else:
        # If entering this branch, set a = b = 0 and choose shoup's algorithm
        a = 0
        b = 0
    return (a,b)

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