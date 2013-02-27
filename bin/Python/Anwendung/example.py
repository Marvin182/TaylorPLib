from TaylorPLib import Matrix, Polynomial
import timeit
import random
import os


# Fill Matrix with random values

random.seed()

def randomPolynomial(order):
	return Polynomial([random.randint(0, 100) / 10 for i in range(order + 1)])

def fillWithRandoms(matrix):
	for i in range(matrix.nrows()):
		for j in range(matrix.ncols()):
			matrix[[i, j]] = randomPolynomial(matrix.order())


# Show operations on Polynomial

def showPolynomial():
	p1 = Polynomial([2, 4, 3])
	p2 = Polynomial([2, 4, 3])

	p3 = p1 + p2
	p4 = p1 * p2

	print("p1:\t\t", p1)
	print("p2:\t\t", p2)
	print("p1 + p2:\t", p3)
	print("p1 * p2:\t", p4)


# Show operations on Matrix

def showMatrix():
	M1 = Matrix(5, 5, 1)
	M1[[0, 0]] = Polynomial([1, 2])
	M1[[0, 1]] = Polynomial([3, 4])
	M1[[0, 2]] = Polynomial([5, 6])
	M1[[0, 3]] = Polynomial([1, 2])
	M1[[0, 4]] = Polynomial([5, 6])

	M1[[1, 0]] = Polynomial([1, 2])
	M1[[1, 1]] = Polynomial([3, 4])
	M1[[1, 2]] = Polynomial([5, 6])
	M1[[1, 3]] = Polynomial([1, 2])
	M1[[1, 4]] = Polynomial([5, 6])

	M1[[2, 0]] = Polynomial([1, 2])
	M1[[2, 1]] = Polynomial([3, 4])
	M1[[2, 2]] = Polynomial([5, 6])
	M1[[2, 3]] = Polynomial([1, 2])
	M1[[2, 4]] = Polynomial([5, 6])

	M1[[3, 0]] = Polynomial([1, 2])
	M1[[3, 1]] = Polynomial([3, 4])
	M1[[3, 2]] = Polynomial([5, 6])
	M1[[3, 3]] = Polynomial([1, 2])
	M1[[3, 4]] = Polynomial([5, 6])

	M1[[4, 0]] = Polynomial([1, 2])
	M1[[4, 1]] = Polynomial([3, 4])
	M1[[4, 2]] = Polynomial([5, 6])
	M1[[4, 3]] = Polynomial([1, 2])
	M1[[4, 4]] = Polynomial([5, 6])

	# lets do it a bit shorter this time
	M2 = Matrix(5, 5, [Polynomial([0, 0]), Polynomial([0, 0]), Polynomial([0, 0]), Polynomial([2, 1]), Polynomial([2, 1]), 
		Polynomial([0, 0]), Polynomial([0, 0]), Polynomial([0, 0]), Polynomial([2, 1]), Polynomial([2, 1]), 
		Polynomial([0, 0]), Polynomial([0, 0]), Polynomial([0, 0]), Polynomial([2, 1]), Polynomial([2, 1]), 
		Polynomial([0, 0]), Polynomial([0, 0]), Polynomial([0, 0]), Polynomial([2, 1]), Polynomial([2, 1]), 
		Polynomial([0, 0]), Polynomial([0, 0]), Polynomial([0, 0]), Polynomial([2, 1]), Polynomial([2, 1])])

	print("M1:", M1)
	print("M2:", M2)

	M3 = M1 + M2
	print("M1 + M2:", M3)

	M4 = M1 * M2
	print("M1 * M2:", M4)

	M5 = M1 * 1.5
	print("M1 * 1.5:", M5)


# Show speed advantage of the otimized matrix multiplication

def setup():
	alpha = random.random()
	beta = random.random()
	m = 10
	n = 10
	order = 5

	A = Matrix(m, n, order)
	B = Matrix(n, n, order)
	C = Matrix(m, n, order)

	fillWithRandoms(A)
	fillWithRandoms(B)
	fillWithRandoms(C)

	return [alpha, beta, A, B, C]


[alpha, beta, A, B, C] = setup()

def normalMatrixMultiplication():
	D = A * B * alpha + C * beta

def optimizedMatrixMultiplication():
	C.mmCaABbC(alpha, beta, A, B)

def showSpecialMatrixMultiplication():
	normal = timeit.timeit("normalMatrixMultiplication()", setup="from __main__ import normalMatrixMultiplication", number=100)
	optimized = timeit.timeit("optimizedMatrixMultiplication()", setup="from __main__ import optimizedMatrixMultiplication", number=100)
	faster = 100 * (1 - optimized / normal)
	print("Die optimierte Matrixmultiplikation war", faster, "% schneller.")


# Run demo

clear = lambda: os.system('cls')
clear()

input("Dies ist ein kleines Programm um die Funktionalität der Klassenbibliothek in Python\nzu demonstrieren. Drücken Sie Enter um zu beginnen.")

clear()
print("\nRechnen mit Polynomen:")
showPolynomial()
input("\nDrücken Sie Enter um fortzufahren.")

clear()
print("\nRechnen mit Matrizen mit Polynomen:")
showMatrix()
input("Drücken Sie Enter um fortzufahren.")

clear()
print("\nOptimierte Matrizenmultiplikation:\n")
showSpecialMatrixMultiplication()

input("\nDas wars. Zum Beenden drücken Sie ein letztes Mal Enter.\n\nVielen Dank für Ihre Aufmerksamkeit.")
