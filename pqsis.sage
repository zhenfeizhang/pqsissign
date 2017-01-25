## system parameters #############

P.<x> = PolynomialRing(ZZ)
N = 443
q = 2^13
F = P(x^N-1)
d = ZZ(round(N/3))
p = 3
p_inv = 1/3 %q
sigma = 90
##################################

def poly_gen(n,d):
  a = P(1)
  i = 0
  while (i<d):
    index = ZZ.random_element(1,n) 
    if (a[index] == 0):
      a = a + P(x^index)
      i = i + 1
  i = 0    
  while (i<d):    
    index = ZZ.random_element(1,n)   
    if (a[index] == 0):
      a = a - P(x^index)
      i = i + 1  
  
  return a
  
def inverse_poly_gen (n,d, q):  
  f = poly_gen(n,d)
  a,b,c = xgcd (f, F)
#  print f,F,a,b,c 
  while (gcd(ZZ(a),q)!=1):
    f = poly_gen(n,d)
    a,b,c = xgcd (f, F)
#    print f,F , a,b,c
    
#  print "?"  
  a_inv = 1/ZZ(a)%q
  f_inv = b*a_inv%q
  
  return f, f_inv    
  
  
## key generation ################
f, f_inv = inverse_poly_gen(N, d, q)
g, g_inv = inverse_poly_gen(N, d, p)
h = g*f_inv*p_inv%F%q

##################################

####up, vp is the digest #########
####of the message ###############

up = (ZZ^N).random_element(x = -1, y = 2)
vp = (ZZ^N).random_element(x = -1, y = 2)
##################################

## signing  ######################
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler

D = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)
r = vector([D() for _ in range (N)])

u1 = p*r+up
v1 = P(u1.list())*h%F%q
a = (P(vp.list())-v1)*g_inv%F%p
v = v1+ a*g%F
b = P(r.list()) + a*f%F
##TODO: Rejection sampling on u,v
##################################


def sign(up, vp, f, g_inv, h):
  D = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)
  r = vector([D() for _ in range (N)])

  u1 = p*r+up
  v1 = P(u1.list())*h%F%q
  a = (P(vp.list())-v1)*g_inv%F%p
  v = v1+ a*g%F
  b = P(r.list()) + a*f%F
  return b


## verify  #######################
u = p*b+P(up.list())
v = u*h%F%q
if( (v - P(vp.list()))%p == 0):
  print "accept"

##################################

