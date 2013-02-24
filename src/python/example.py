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



print("Polynome: \n")
# showPolynom()
# Console.ReadKey()
# Console.Clear()

print("Matrizen mit Polynomen: \n")
# showMatrix()
# Console.ReadKey()
# Console.Clear()

print("Spezielle Matrizenmultiplikation: \n")
# showSpecialMatrixMultiplication()
# Console.ReadKey()



def showPolynom():
	p1 = Polynomial([2, 4, 3])
	p2 = Polynomial([2, 4, 3])

	p3 = p1 + p2
	p4 = p1 * p2
	p5 = p4 / p2

	print("p1:", p1)
	print("p2:", p2)
	print("p1 + p2:", p3)
	print("p1 * p2:", p4)
	print("p1 * p2 / p2:", p5)

# def showMatrix():
# 	Polynomial[,] Pa = new Polynomial[5, 5]
# 	Pa[0, 0] = Polynomial([1, 2])
# 	Pa[0, 1] = Polynomial([3, 4])
# 	Pa[0, 2] = Polynomial([5, 6])
# 	Pa[0, 3] = Polynomial([1, 2])
# 	Pa[0, 4] = Polynomial([5, 6])

# 	Pa[1, 0] = Polynomial([1, 2])
# 	Pa[1, 1] = Polynomial([3, 4])
# 	Pa[1, 2] = Polynomial([5, 6])
# 	Pa[1, 3] = Polynomial([1, 2])
# 	Pa[1, 4] = Polynomial([5, 6])

# 	Pa[2, 0] = Polynomial([1, 2])
# 	Pa[2, 1] = Polynomial([3, 4])
# 	Pa[2, 2] = Polynomial([5, 6])
# 	Pa[2, 3] = Polynomial([1, 2])
# 	Pa[2, 4] = Polynomial([5, 6])

# 	Pa[3, 0] = Polynomial([1, 2])
# 	Pa[3, 1] = Polynomial([3, 4])
# 	Pa[3, 2] = Polynomial([5, 6])
# 	Pa[3, 3] = Polynomial([1, 2])
# 	Pa[3, 4] = Polynomial([5, 6])

# 	Pa[4, 0] = Polynomial([1, 2])
# 	Pa[4, 1] = Polynomial([3, 4])
# 	Pa[4, 2] = Polynomial([5, 6])
# 	Pa[4, 3] = Polynomial([1, 2])
# 	Pa[4, 4] = Polynomial([5, 6])

# 	Polynomial[,] Pb = new Polynomial[5, 5]
# 	Pb[0, 0] = Polynomial([0, 0])
# 	Pb[0, 1] = Polynomial([0, 0])
# 	Pb[0, 2] = Polynomial([0, 0])
# 	Pb[0, 3] = Polynomial([2, 1])
# 	Pb[0, 4] = Polynomial([2, 1])

# 	Pb[1, 0] = Polynomial([0, 0])
# 	Pb[1, 1] = Polynomial([0, 0])
# 	Pb[1, 2] = Polynomial([0, 0])
# 	Pb[1, 3] = Polynomial([2, 1])
# 	Pb[1, 4] = Polynomial([2, 1])

# 	Pb[2, 0] = Polynomial([0, 0])
# 	Pb[2, 1] = Polynomial([0, 0])
# 	Pb[2, 2] = Polynomial([0, 0])
# 	Pb[2, 3] = Polynomial([2, 1])
# 	Pb[2, 4] = Polynomial([2, 1])

# 	Pb[3, 0] = Polynomial([0, 0])
# 	Pb[3, 1] = Polynomial([0, 0])
# 	Pb[3, 2] = Polynomial([0, 0])
# 	Pb[3, 3] = Polynomial([2, 1])
# 	Pb[3, 4] = Polynomial([2, 1])

# 	Pb[4, 0] = Polynomial([0, 0])
# 	Pb[4, 1] = Polynomial([0, 0])
# 	Pb[4, 2] = Polynomial([0, 0])
# 	Pb[4, 3] = Polynomial([2, 1])
# 	Pb[4, 4] = Polynomial([2, 1])


# 	Matrix M1 = new Matrix(5, 5, Pa)
# 	Matrix M2 = new Matrix(5, 5, Pb)
# 	Matrix M3 = M1 + M2
# 	Matrix M4 = M1 * M2

# 	print("M1: \n" + M1.ToString())
# 	print("M2: \n" + M2.ToString())
# 	print("M1 + M2: \n" + M3.ToString())
# 	print("M1 * M2: \n" + M4.ToString())

# def showSpecialMatrixMultiplication():
# 	Polynomial[,] Pa = new Polynomial[3, 3]
# 	Pa[0, 0] = Polynomial([1, 2])
# 	Pa[0, 1] = Polynomial([3, 4])
# 	Pa[0, 2] = Polynomial([5, 6])

# 	Pa[1, 0] = Polynomial([7, 8])
# 	Pa[1, 1] = Polynomial([9, 0])
# 	Pa[1, 2] = Polynomial([1, 2])

# 	Pa[2, 0] = Polynomial([3, 4])
# 	Pa[2, 1] = Polynomial([5, 6])
# 	Pa[2, 2] = Polynomial([7, 8])

# 	Polynomial[,] Pb = new Polynomial[3, 3]
# 	Pb[0, 0] = Polynomial([3, 4])
# 	Pb[0, 1] = Polynomial([5, 6])
# 	Pb[0, 2] = Polynomial([7, 8])

# 	Pb[1, 0] = Polynomial([1, 2])
# 	Pb[1, 1] = Polynomial([3, 4])
# 	Pb[1, 2] = Polynomial([5, 6])

# 	Pb[2, 0] = Polynomial([7, 8])
# 	Pb[2, 1] = Polynomial([9, 0])
# 	Pb[2, 2] = Polynomial([1, 2])

# 	Polynomial[,] Pc = new Polynomial[3, 3]
# 	Pc[0, 0] = Polynomial([7, 8])
# 	Pc[0, 1] = Polynomial([9, 0])
# 	Pc[0, 2] = Polynomial([1, 2])

# 	Pc[1, 0] = Polynomial([3, 4])
# 	Pc[1, 1] = Polynomial([5, 6])
# 	Pc[1, 2] = Polynomial([7, 8])

# 	Pc[2, 0] = Polynomial([1, 2])
# 	Pc[2, 1] = Polynomial([3, 4])
# 	Pc[2, 2] = Polynomial([5, 6])


# 	Matrix M1 = new Matrix(3, 3, Pa)
# 	Matrix M2 = new Matrix(3, 3, Pb)
# 	Matrix M3 = new Matrix(3, 3, Pc)

# 	print("( M1 * M2 * alpha ) + ( M3 * beta ) where: ")
# 	print("M1: \n" + M1.ToString())
# 	print("M2: \n" + M2.ToString())
# 	print("M3: \n" + M3.ToString())
# 	print("alpha: 2")
# 	print("beta: 2\n")

# 	DateTime start = DateTime.Now
# 	String temp = ((M1 * M2 * 2) + (M3 * 2)).ToString()
# 	print(temp)
# 	DateTime end = DateTime.Now
# 	print("Time needed: " + (end - start).ToString())
# 	print("This was normal mode...")
# 	start = DateTime.Now
# 	M3.mmCaABbC(2, 2, M1, M2)
# 	String temp2 = M3.ToString()
# 	print(temp2)
# 	end = DateTime.Now
# 	print("Time needed: " + (end - start).ToString())
# 	print("This was method mmCaABbC...")

#     if (temp.Equals(temp2)):
# 		print("Both Are Equal...")
#     else:
# 	    print("Both Are NOT Equal. Something went wrong...")
