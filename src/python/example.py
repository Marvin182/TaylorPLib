from TaylorPLib import Matrix, Polynomial
from array import array

p1 = Polynomial([2.0, 6.0, 2.0, 4.0])
print("p1=", p1)
print("p1[1] = 1.0")
p1[1] = 1.0
print("p1=", p1)

p2 = Polynomial(1)
p2.setCoeffs([3.0, 5.0])
print("p2=", p2)


M = Matrix(2, 3)
print("M=", M)
print("M(1, 1) = 2 + 3x")
M[[1, 1]] = Polynomial([2.0, 3.0])
print(M)
