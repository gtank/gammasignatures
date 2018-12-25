import binascii
import hashlib

# Twisted Edwards curve -1*x^2+y^2=1+d*x^2*y^2

# ed25519 parameters
K = GF(2^255 - 19)
a = K(-1)
d = K(-121665/121666)
base_y = K(4/5)
base_x = ((1 - base_y^2) / (-1 - d*(base_y^2))).sqrt()

# Convert Edwards to Montgomery
A = 2*(a+d)/(a-d)
B = 4/(a-d)
base_u = (1+base_y)/(1-base_y)
base_v = base_u/base_x

# Convert Montgomery to Weierstrass
a = (3-A^2)/(3*B^2)
b = (2*A^3-9*A)/(27*B^3)
x0, y0 = (base_u+A/3)/B, base_v/B

# Short Weierstrass form:
# y^2 = x^3 + 42204101795669822316448953119945047945709099015225996174933988943478124189485*x +
# 13148341720542919587570920744190446479425344491440436116213316435534172959396

E = EllipticCurve(GF(2^255-19), [a, b])
Scalars = GF(E.cardinality() / 8)
P = E(x0, y0)

# gamma signatures
def KeyGen():
	x = Scalars.random_element()
	X = Integer(x) * P
	return (x,X)

def H_d(point):
	return Scalars(point[0])

def H_e(X, m):
	h = hashlib.sha256()
	h.update(X.dumps())
	h.update(m)
	return Scalars(Integer(binascii.hexlify(h.digest()), base=16))

def Sign(X, x, m):
	r = Scalars.random_element()
	A = Integer(r)*P
	d = H_d(A)
	e = H_e(X, m)
	z = r*d - e*x
	return (d, z)

def Verify(X, m, d, z):
    e = H_e(X, m)
    d_inv = inverse_mod(Integer(d), Scalars.cardinality())
    A = Integer((z*d_inv))*P+Integer((e*d_inv))*X
    if d == H_d(A):
        return True
    else:
        return False

x, X = KeyGen()
d, z = Sign(X, x, "Hello, world")
print(Verify(X, "Hello, world", d, z))
