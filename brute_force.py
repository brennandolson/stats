import numpy as np

cache = {}
def C(n,k):
  '''n Choose k'''
  if (n,k) in cache:
    return cache[(n,k)]
  elif(n, n - k) in cache:
    return cache[(n, n - k)]
  else:
    ret = 1
    for i in range(1, k + 1):
      ret *= (n + 1 - i)
      ret //= i
    cache[(n, k)] = ret
    return ret

def half_paths(n, k):
  '''generates the valid paths of
     R and D which walk halfway through
     an (n+1)-by-(n+1) grid starting
     in the upper right corner
     without crossing the diagonal'''
  if n == 0:
    yield ''
    return
  
  if k > 0:
    for path in half_paths(n - 1, k - 1):
      yield 'D' + path
  
  for path in half_paths(n - 1, k + 1):
    yield 'R' + path
  
  return

def full_paths(n):
  '''generates the valid SYMMETRIC paths of
     R and D which walk through
     an (n+1)-by-(n+1) grid starting
     in the upper right corner
     without crossing the diagonal'''
  for path in half_paths(n, 0):
    suffix = ''
    for l in path:
      if l == 'R':
        suffix = 'L' + suffix
      else:
        suffix = 'R' + suffix
    yield path + suffix
  return

def get_polynomials_from_path(path, n):
  '''Given a path through an n-by-n grid,
     interprets the walk as a sequence of polynomials
     in the context of the CI problem'''
  polys = []
  basic_polys = []
  p = np.poly1d([1,0])   # y = x
  q = np.poly1d([-1, 1]) # y = 1 - x
  for i in range(n):
    basic_polys.append(C(n - 1, i) * p**i * q**(n - 1 - i))
  
  i, j = 0, 0
  for c in path:
    polys.append(sum(basic_polys[i:j + 1]))
    if c == 'R':
      j += 1
    else:
      i += 1
  polys.append(basic_polys[-1])
  return polys 

def newton(p, epsilon = 0.0001):
  der = np.polyder(p)
  x0 = 0.5
  x1 = 0.5 - np.polyval(p, 0.5) / np.polyval(der, 0.5)
  while abs(x1 - x0) > epsilon:
    x0 = x1
    x1 = x0 - np.polyval(p, x0) / np.polyval(der, x0)
  return x1

polys = get_polynomials_from_path("RRDRDD", 4)

def jump_point(f, g, a):
  '''where to jump from f(x) to g(x) to minimize
     difference of (h(x) - a)^2 where h takes
     the value of either f or g depending on 
     the jump point computed here'''
  diff_of_squares = (f - a)**2 - (g - a)**2
  return newton(diff_of_squares)

for i in range(4):
  print (jump_point(polys[i], polys[i + 1], 0.95))
