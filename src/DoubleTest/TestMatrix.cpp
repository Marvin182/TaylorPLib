#include <stdio.h>
#include <gtest/gtest.h>
#include "TaylorPLib.h"

using namespace std;
using namespace LibMatrix;

#define T(n) Polynomial() * (n)

namespace LibMatrix {
	// needed by GTest to print matrices if some ASSERT_EQ(matrix1, matrix2) failed
	std::ostream& operator<<(std::ostream &out, const Matrix &m)
	{
		out << setiosflags(ios::fixed) << setprecision(2);
		for (int i = 0; i < m.nrows(); i++)
		{
			out << '\n';
			for (int j = 0; j < m.ncols(); j++)
			{
				out << m.get(i, j) << '\t';
			}
		}
		return out; 
	}

	bool operator==(double x, const Polynomial &p)
	{
		return p.isConst() && p[0] == x;
	}
}

/*
 * FIXTURES
 */
class MatrixOperator: public ::testing::Test
{
	protected:
  		Matrix A, B, AplusB, minusA, AminusB, twoA, AB, BA;

		MatrixOperator()
		{
			// A
			Polynomial a[] = {
				T(1), T(2),
				T(3), T(4)
			};
			A = Matrix(2, 2, a);

			// B
			Polynomial b[] = {
				T(5), T(2),
				T(0), T(3)
			};
			B = Matrix(2, 2, b);

			// A + B
			Polynomial aplusb[] = {
				T(6), T(4),
				T(3), T(7) 
			};
			AplusB = Matrix(2, 2, aplusb);

			// - A
			Polynomial minusa[] = {
				T(-1), T(-2),
				T(-3), T(-4)
			};
			minusA = Matrix(2, 2, minusa);

			// A - B
			Polynomial aminusb[] = {
				T(-4), T(0),
				 T(3), T(1)
			};
			AminusB = Matrix(2, 2, aminusb);

			// T(2) * A
			Polynomial twoa[] = {
				T(2), T(4),
				T(6), T(8)
			};
			twoA = Matrix(2, 2, twoa);

			// A * B
			Polynomial ab[] = {
				 T(5),  T(8),
				T(15), T(18)
			};
			AB = Matrix(2, 2, ab);			
			
			// B * A
			Polynomial ba[] = {
				T(11), T(18),
				 T(9), T(12)
			};
			BA = Matrix(2, 2, ba);
		}
};

class MatrixMultiplication: public ::testing::Test
{
	protected:
		Matrix A, B, C, I;
		double alpha, beta;

		Polynomial getRandomNumber() { return T(rand() % 10); }
		void fillWithRandoms(Matrix &m)
		{
			for (int i = 0; i < m.nrows(); i++)
			{
				for (int j = 0; j < m.ncols(); j++)
				{
					m(i, j) = getRandomNumber();
				}
			}
		}

		MatrixMultiplication():
			alpha(2.0),
			beta(3.0)
		{
			A = Matrix(5, 5);
			fillWithRandoms(A);

			B = Matrix(5, 5);
			fillWithRandoms(B);

			C = Matrix(5, 5);
			fillWithRandoms(C);

			I = Matrix(5, 5);
			I.set2Id();
		}
};

class MatrixMethods: public ::testing::Test
{
	protected:
		Matrix I, Z, A, B, C, CZ, CId, C3;

		MatrixMethods()
		{
			Polynomial i[] = {
				T(1), T(0), T(0),
				T(0), T(1), T(0),
				T(0), T(0), T(1)
			};
			I = Matrix(3, 3, i);

			Polynomial z[] = {
				T(0), T(0), T(0),
				T(0), T(0), T(0),
				T(0), T(0), T(0)
			};
			Z = Matrix(3, 3, z);

			Polynomial a[] = {
				T(1), T(2), T(3),
				T(4), T(5), T(6),
				T(7), T(8), T(9)
			};
			A = Matrix(3, 3, a);

			Polynomial b[] = {
				T(1),  T(2),  T(3),  T(4),
				T(5),  T(6),  T(7),  T(8),
				T(9), T(10), T(11), T(12)
			};
			B = Matrix(3, 4, b);

			Polynomial c[] = {
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5)
			};
			C = Matrix(8, 7, c);

			Polynomial cz[] = {
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(5), T(0), T(0), T(5), T(5), T(5),
				T(5), T(5), T(0), T(0), T(5), T(5), T(5),
				T(5), T(5), T(0), T(0), T(5), T(5), T(5),
				T(5), T(5), T(0), T(0), T(5), T(5), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5)
			};
			CZ = Matrix(8, 7, cz);

			Polynomial cid[] = {
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(5), T(1), T(0), T(0), T(0), T(5),
				T(5), T(5), T(0), T(1), T(0), T(0), T(5),
				T(5), T(5), T(0), T(0), T(1), T(0), T(5),
				T(5), T(5), T(0), T(0), T(0), T(1), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5)
			};
			CId = Matrix(8, 7, cid);

			Polynomial c3[] = {
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(3), T(3), T(3), T(3), T(5), T(5),
				T(5), T(3), T(3), T(3), T(3), T(5), T(5),
				T(5), T(3), T(3), T(3), T(3), T(5), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5),
				T(5), T(5), T(5), T(5), T(5), T(5), T(5)
			};
			C3 = Matrix(8, 7, c3);
		}
};

class MatrixExceptions: public ::testing::Test
{
protected:
	Matrix A,B;

	MatrixExceptions()
	{
		Polynomial a[] = {
			T(1), T(6), T(0),
			T(2), T(7), T(0),
			T(3), T(8), T(0)
		};
		Polynomial b[] = {
			T(1), T(6),
			T(2), T(7),
			T(3), T(8)
		};

		A = Matrix(3, 3, a);
		B = Matrix(3, 2, b);
	}
};

/*
 * CONSTRUCTOR TESTS
 */
TEST(MatrixConstructor, default_constructor)
{
	// default constructor should create a zero initialized T(1)x1 matrix 
	Matrix d;
	ASSERT_EQ(1, d.nrows());
	ASSERT_EQ(1, d.ncols());
	ASSERT_EQ(T(0), d(0, 0));
}

TEST(MatrixConstructor, regular_constructor)
{	
	// regular constructor
	Matrix r(2, 3);
	ASSERT_EQ(2, r.nrows());
	ASSERT_EQ(3, r.ncols());
	ASSERT_EQ(0, r(0, 0));
	ASSERT_EQ(0, r(1, 2));
}

TEST(MatrixConstructor, test_and_copy_constructor)
{
	// test constructor
	Polynomial values[] = {
		T(1), T(2),	// first row
		T(3), T(4),	// second row
		T(5), T(6)	// third row
	};
	Matrix t(3, 2, values); // rows, columns, values
	// t.print("t");
	ASSERT_EQ(3, t.nrows());
	ASSERT_EQ(2, t.ncols());
	ASSERT_EQ(1, t(0, 0));
	ASSERT_EQ(2, t(0, 1));
	ASSERT_EQ(3, t(1, 0));
	ASSERT_EQ(4, t(1, 1));
	ASSERT_EQ(5, t(2, 0));
	ASSERT_EQ(6, t(2, 1));

	// copy constructor
	Matrix c(t);
	// c should be a deep copy of t
	t(0, 1) = T(8);
	t(2, 0) = T(9);
	ASSERT_EQ(8, t(0, 1));
	ASSERT_EQ(9, t(2, 0));

	ASSERT_EQ(3, c.nrows());
	ASSERT_EQ(2, c.ncols());
	ASSERT_EQ(1, c(0, 0));
	ASSERT_EQ(2, c(0, 1));
	ASSERT_EQ(3, c(1, 0));
	ASSERT_EQ(4, c(1, 1));
	ASSERT_EQ(5, c(2, 0));
	ASSERT_EQ(6, c(2, 1));
}

TEST_F(MatrixOperator, assignment)
{
	B = A;

	ASSERT_EQ(2, B.nrows());
	ASSERT_EQ(2, B.ncols());

	for (int i = 0; i < A.nrows(); i++)
	{
		for (int j = 0; j < A.ncols(); j++)
		{
			ASSERT_EQ(A(i, j), B(i, j));
		}
	}
}

TEST_F(MatrixOperator, comparsion)
{
	ASSERT_TRUE(A == A);	
	ASSERT_FALSE(A == B);

	ASSERT_FALSE(A != A);
	ASSERT_TRUE(A != B);

	Matrix asA(A);
	ASSERT_TRUE(A == asA);

	Matrix biggerA(2, 3);
	for (int i = 0; i < A.nrows(); i++)
	{
		for (int j = 0; j < A.ncols(); j++)
		{
			biggerA(i, j) = A(i, j);
		}
	}
	ASSERT_FALSE(A == biggerA);

	// ASSERT_EQ and ASSERT_EQ should use the comparision operators
	ASSERT_EQ(A, asA);
	ASSERT_NE(A, B);
}

TEST_F(MatrixOperator, plus)
{
	ASSERT_EQ(AplusB, A + B);

	A += B;
	ASSERT_EQ(AplusB, A);
}

TEST_F(MatrixOperator, unary)
{
	ASSERT_EQ(minusA, -A);
}

TEST_F(MatrixOperator, minus)
{
	ASSERT_EQ(AminusB, A - B);

	A -= B;
	ASSERT_EQ(AminusB, A);
}

TEST_F(MatrixOperator, timesScalar)
{
	ASSERT_EQ(twoA, A * 2.0);

	A *= 2.0;
	ASSERT_EQ(twoA, A);
}

TEST_F(MatrixOperator, timesMatrix)
{
	ASSERT_EQ(AB, A * B);
	ASSERT_EQ(BA, B * A);
}

/*
 * SPECIAL MATRIX MULTIPLICATIONS
 */
TEST_F(MatrixMultiplication, mmCaABbC)
{
	Matrix expect = (A * B * alpha) + (C * beta);
	C.mmCaABbC(alpha, beta, A, B);

	ASSERT_EQ(expect, C);
}

TEST_F(MatrixMultiplication, bmmCaABbC)
{
   	// special form for A
	Polynomial b[] = {
		T(1), T(2), T(3), T(0), T(0),
		T(4), T(5), T(6), T(0), T(0),
		T(7), T(8), T(9), T(0), T(0),
		T(0), T(0), T(0), T(1), T(0),
		T(0), T(0), T(0), T(0), T(1)
	};

	B = Matrix(5, 5, b);

	Matrix expect = (A * B * alpha) + (C * beta);
	C.bmmCaABbC(3, 3, alpha, beta, A, B);

	ASSERT_EQ(expect, C);	
}

TEST_F(MatrixMultiplication, mmCasABbC)
{
   	// special form for A
	Polynomial a[] = {
		T(0), T(0), T(0), T(0), T(0),
		T(0), T(0), T(0), T(0), T(0),
		T(0), T(0), T(0), T(0), T(0),
		T(1), T(2), T(3), T(4), T(5),
		T(5), T(4), T(3), T(2), T(1)
	};
	A = Matrix(5, 5, a);

	Matrix expect = (A * B * alpha) + (C * beta);
	C.mmCasABbC(2, alpha, beta, A, B);
	ASSERT_EQ(expect, C);	
}

TEST_F(MatrixMultiplication, mmCaAsBbC)
{
   	// special form for B
	Polynomial b[] = {
		T(0), T(0), T(0), T(1), T(5),
		T(0), T(0), T(0), T(2), T(4),
		T(0), T(0), T(0), T(3), T(3),
		T(0), T(0), T(0), T(4), T(2),
		T(0), T(0), T(0), T(5), T(1)
	};
	B = Matrix(5, 5, b);

	Matrix expect = (A * B * alpha) + (C * beta);
	C.mmCaAsBbC(2, alpha, beta, A, B);
	ASSERT_EQ(expect, C);	
}

TEST_F(MatrixMultiplication, mmCaAUTBPbC)
{
	Polynomial bUT[] = {
		T(1), T(2), T(3), T(4), T(5),
		T(0), T(6), T(7), T(8), T(9),
		T(0), T(0), T(1), T(2), T(3),
		T(0), T(0), T(0), T(4), T(5),
		T(0), T(0), T(0), T(0), T(6)
	};
	Matrix B_UT = Matrix(5, 5, bUT);
	
	// special form for B
	Polynomial b[] = {
		T(1), T(3), T(2), T(4), T(5),
		T(0), T(7), T(6), T(8), T(9),
		T(0), T(1), T(0), T(2), T(3),
		T(0), T(0), T(0), T(4), T(5),
		T(0), T(0), T(0), T(0), T(6)
	};
	B = Matrix(5, 5, b);
	int piv[] = { 0, 2, 1, 3, 4 } ;

	Matrix expect = (A * B_UT * alpha) + (C * beta);
	C.mmCaAUTBPbC(alpha, beta, A, B, (int*) piv);
	ASSERT_EQ(expect, C);
}

TEST_F(MatrixMultiplication, mmCaAATbC)
{
	Matrix at = A.asTranspose();
    Matrix expect = (A * at * alpha) + (C * beta);
	C.mmCaAATbC(alpha, beta, A);
	ASSERT_EQ(expect, C);
}

TEST_F(MatrixMultiplication, mmCaATAbC)
{
	Matrix at = A.asTranspose();
    Matrix expect = (at * A * alpha) + (C * beta);
	C.mmCaATAbC(alpha, beta, A);
	ASSERT_EQ(expect, C);
}

TEST_F(MatrixMultiplication, mmCaATBbC)
{
	Matrix at = A.asTranspose();
    Matrix expect = (at * B * alpha) + (C * beta);
	C.mmCaATBbC(alpha, beta, A, B);
	ASSERT_EQ(expect, C);
}

TEST_F(MatrixMultiplication, mmCaATBPbC)
{
	Polynomial a1[] = {
		T(1), T(2), T(3), T(4), T(5),
		T(6), T(7), T(8), T(9), T(1),
		T(0), T(0), T(0), T(0), T(0),
		T(0), T(0), T(0), T(0), T(0),
		T(0), T(0), T(0), T(0), T(0)
	};
	
	Polynomial b1[] = {
		T(1), T(2), T(3), T(4), T(5),
		T(6), T(7), T(8), T(9), T(1),
		T(0), T(0), T(0), T(0), T(0),
		T(0), T(0), T(0), T(0), T(0),
		T(0), T(0), T(0), T(0), T(0)
	};
	
	Polynomial b2[] = {
		T(1), T(2), T(3), T(5), T(4),
		T(6), T(7), T(8), T(1), T(9),
		T(0), T(0), T(0), T(0), T(0),
		T(0), T(0), T(0), T(0), T(0),
		T(0), T(0), T(0), T(0), T(0)
	};

	A = Matrix(5, 5, a1);
	B = Matrix(5, 5, b1);
	Matrix at = A.asTranspose();
	Matrix expect = (at * B * alpha) + (C * beta);

	int piv[] = { 0, 1, 2, 4, 3 } ;
	int *pointer = piv;

	B = Matrix(5, 5, b2);
	C.mmCaATBPbC(alpha, beta, A, B, pointer);
	ASSERT_EQ(expect, C);
}

TEST_F(MatrixMultiplication, mmCaABTbC1)
{
	Matrix bt = B.asTranspose();
    Matrix expect = (A * bt * alpha) + (C * beta);
	C.mmCaABTbC(alpha, beta, A, B);
	ASSERT_EQ(expect, C);
}

TEST_F(MatrixMultiplication, mmCaABTbC_up)
{
	Polynomial b[] = {
		T(1), T(6), T(0), T(0), T(0),
		T(2), T(7), T(0), T(0), T(0),
		T(3), T(8), T(0), T(0), T(0),
		T(4), T(9), T(0), T(0), T(0),
		T(5), T(1), T(0), T(0), T(0)
	};

	B = Matrix(5, 5, b);
	Matrix bt = B.asTranspose();
    Matrix expect = (A * bt * alpha) + (C * beta);

	C.mmCaABTbC(2, true, alpha, beta, A, B);
	ASSERT_EQ(expect, C);
}

TEST_F(MatrixMultiplication, mmCaABTbC_down)
{
	Polynomial b[] = {
		T(0), T(0), T(0), T(1), T(6),
		T(0), T(0), T(0), T(2), T(7),
		T(0), T(0), T(0), T(3), T(8),
		T(0), T(0), T(0), T(4), T(9),
		T(0), T(0), T(0), T(5), T(1)
	};

	B = Matrix(5, 5, b);
	Matrix bt = B.asTranspose();
    Matrix expect = (A * bt * alpha) + (C * beta);

	C.mmCaABTbC(2, false, alpha, beta, A, B);
	ASSERT_EQ(expect, C);
}

TEST_F(MatrixMultiplication, bmmCaABTbC)
{
	// case T(1)
	Polynomial a[] = {
		T(1), T(2), T(3), T(0), T(0),
		T(4), T(5), T(6), T(0), T(0),
		T(7), T(8), T(9), T(0), T(0),
		T(0), T(0), T(0), T(1), T(0),
		T(0), T(0), T(0), T(0), T(1)
	};
	A = Matrix(5, 5, a);
	
	Matrix bt = B.asTranspose();
	Matrix expect = (A * bt * alpha) + (C * beta);
	Matrix tempC = C;

	C.bmmCaABTbC(3, 3, alpha, beta, A, B);
	ASSERT_EQ(expect, C);

	tempC.bmmCaABTbC(2, 2, alpha, beta, A, B);
	ASSERT_NE(expect, tempC);

	// case T(2)
	Polynomial a2[] = {
		T(1), T(2), T(3), T(0),
		T(4), T(5), T(6), T(0),
		T(0), T(0), T(0), T(1)
	};
	Polynomial b2[] = {
		T(1), T(2), T(3),
		T(4), T(5), T(6),
		T(0), T(0), T(0),
		T(0), T(0), T(0)
	};

	A = Matrix(3, 4, a2);
	B = Matrix(4, 3, b2);
	C = Matrix(3, 3);
	// Exception, da B transposed wird und da die Mulitplikation (3 x 4) * (3 x 4) nicht möglich ist
	ASSERT_THROW(C.bmmCaABTbC(2, 3, alpha, beta, A, B), MathException);
}

TEST_F(MatrixMultiplication, mmCaIBbC)
{
	Matrix expect = (I * B * alpha) + (C * beta);

	C.mmCaIBbC(alpha, beta, B);
	ASSERT_EQ(expect, C);
}
// void Matrix::mmCaIBbC(double alpha, double beta, int *piv, bool rows, const Matrix&B);

TEST_F(MatrixMultiplication, mmCaAIbC)
{
	Matrix expect = (A * I * alpha) + (C * beta);

	C.mmCaAIbC(alpha, beta, A);
	ASSERT_EQ(expect, C);
}
// void Matrix::mmCaAIbC(double alpha, double beta, const Matrix &A, int *piv, bool rows);

/*
 * OTHER MATRIX METHODS
 */
// double Matrix::fnorm();
// int Matrix::tcfnorm(int nrTC, double *fn);
// int Matrix::colnorm(double *c);
// int Matrix::colnormdown(int pos, int *piv, double *c);

TEST_F(MatrixMethods, cpermutem)
{
	int newOrder[] = {3, 0, 2, 1};
	Polynomial b[] = {
		 T(4),  T(1),  T(3),  T(2),
		 T(8),  T(5),  T(7),  T(6),
		T(12),  T(9),  T(11), T(10)
	};
	Matrix newB = Matrix(3, 4, b);

	B.cpermutem((int*) newOrder);
	ASSERT_EQ(newB, B);	
}

TEST_F(MatrixMethods, cpermutem_trans)
{
	int newOrder[] = {1, 3, 2, 0};
	// Beim Transpose wird aus dem Pivot Vektor -> { 3, 0, 2, 1 }

	Polynomial b[] = {
		 T(4),  T(1),  T(3),  T(2),
		 T(8),  T(5),  T(7),  T(6),
		T(12),  T(9),  T(11), T(10)
	};
	Matrix newB = Matrix(3, 4, b);

	B.cpermutem((int*) newOrder, true);
	ASSERT_EQ(newB, B);	
}

TEST_F(MatrixMethods, rpermutem)
{
	int newOrder[] = {1, 2, 0};
	Polynomial b[] = {
		T(5),  T(6),  T(7),  T(8),
		T(9), T(10), T(11), T(12),
		T(1),  T(2),  T(3),  T(4)
	};
	Matrix newB = Matrix(3, 4, b);

	B.rpermutem((int*) newOrder);
	ASSERT_EQ(newB, B);
}

TEST_F(MatrixMethods, transpose)
{
	// test in place transpose
	Polynomial atrans[] = {
		T(1), T(4), T(7),
		T(2), T(5), T(8),
		T(3), T(6), T(9)
	};
	Matrix Atranspose(3, 3, atrans);
	A.transpose();
	ASSERT_EQ(Atranspose, A);

	// test normal transpose
	Polynomial btrans[] = {
		T(1),  T(5),  T(9),
		T(2),  T(6), T(10),
		T(3),  T(7), T(11),
		T(4),  T(8), T(12)
	};
	Matrix Btranspose(4, 3, btrans);
	B.transpose();
	ASSERT_EQ(Btranspose, B);	
}

TEST_F(MatrixMethods, asTranspose)
{
	Polynomial atrans[] = {
		T(1), T(4), T(7),
		T(2), T(5), T(8),
		T(3), T(6), T(9)
	};
	Matrix Atranspose(3, 3, atrans);
	ASSERT_EQ(Atranspose, A.asTranspose());

	Polynomial btrans[] = {
		T(1),  T(5),  T(9),
		T(2),  T(6), T(10),
		T(3),  T(7), T(11),
		T(4),  T(8), T(12)
	};
	Matrix Btranspose(4, 3, btrans);
	ASSERT_EQ(Btranspose, B.asTranspose());
}

TEST_F(MatrixMethods, shift)
{
	// shift would be senseless for double matrices
}

TEST_F(MatrixMethods, isId)
{
	ASSERT_TRUE(I.isId());
	ASSERT_FALSE(Z.isId());
	ASSERT_FALSE(A.isId());
	ASSERT_FALSE(B.isId());		
}
// bool Matrix::isId(double eps);
// bool Matrix::isId(int m1, int m2, int n1, int n2, double eps);

TEST_F(MatrixMethods, isZero)
{
 	ASSERT_FALSE(I.isZero());
	ASSERT_TRUE(Z.isZero());
	ASSERT_FALSE(A.isZero());
	ASSERT_FALSE(B.isZero());
}
// bool Matrix::isZero(double eps);

TEST_F(MatrixMethods, set2Id)
{
	A.set2Id();
	ASSERT_EQ(I, A);
}

TEST_F(MatrixMethods, set2Id2)
{
	C.set2Id(1, 3, 2, 1);
	ASSERT_EQ(CId, C);
}

TEST_F(MatrixMethods, set2IdFromIndices)
{
	C.set2IdFromIndices(1, 4, 2, 5);
	ASSERT_EQ(CId, C);
}

TEST_F(MatrixMethods, set2Zero)
{
	A.set2Zero();
	ASSERT_EQ(Z, A);

	Polynomial zeros[] = {
		T(0), T(0), T(0), T(0),
		T(0), T(0), T(0), T(0),
		T(0), T(0), T(0), T(0)
	};
	Matrix Zeros(3, 4, zeros);
	
	B.set2Zero();
	ASSERT_EQ(Zeros, B);
}

TEST_F(MatrixMethods, set2Zero2)
{
	C.set2Zero(1, 3, 2, 3);
	ASSERT_EQ(CZ, C);
}

TEST_F(MatrixMethods, set2ZeroFromIndices)
{
	C.set2ZeroFromIndices(1, 4, 2, 3);
	ASSERT_EQ(CZ, C);
}

TEST_F(MatrixMethods, set2Val)
{
	Polynomial b42[] = {
		T(42), T(42), T(42), T(42),
		T(42), T(42), T(42), T(42),
		T(42), T(42), T(42), T(42),
	};
	Matrix B42(3, 4, b42);

	B.set2Val(42.0);
	ASSERT_EQ(B42, B);
}

TEST_F(MatrixMethods, set2Val2)
{
	// C.set2Val(2, 3, T(1), T(2), T(3).0);
	ASSERT_EQ(C3, C);
}

TEST_F(MatrixMethods, set2ValFromIndices)
{
	// C.set2ValFromIndices(2, 4, T(1), T(4), T(3).0);
	ASSERT_EQ(C3, C);
}

// int Matrix::trinvm(Matrix &Am);
// int Matrix::trinvm(int r, Matrix &Am);
// bool Matrix::mcompare(Matrix &B);
// bool Matrix::mcompare(Matrix &B, double eps);
// bool Matrix::mcompare(Matrix &B, int r, double eps);

TEST_F(MatrixExceptions, OperatorExceptions)
{
	ASSERT_THROW(A.get(4,3), MathException);
	ASSERT_THROW(A(4,3), MathException);
	ASSERT_THROW(A + B, MathException);
	ASSERT_THROW(A += B, MathException);
	ASSERT_THROW(A - B, MathException);
	ASSERT_THROW(A -= B, MathException);
	ASSERT_THROW(B * A, MathException);
}

TEST_F(MatrixExceptions, FunctionExceptions)
{
	// mmCaABbC
	// Wrong Dimension of receiving Matrix
	Matrix C(3, 3);
	ASSERT_THROW(C.mmCaABbC(1, 1, A, B),MathException);
	C = Matrix(3, 2);
	C.mmCaABbC(1, 1, A, B); // should be okay
	// Wrong Multiplication Dimensions
	ASSERT_THROW(C.mmCaABbC(1, 1, B, A),MathException);

	// bmmCaABbC
	// rows > rows of B
	ASSERT_THROW(C.bmmCaABbC(5, 2, 1, 1, A, B), MathException);
	// cols > cols of B
	ASSERT_THROW(C.bmmCaABbC(2, 5, 1, 1, A, B), MathException);
	// A.cols != B.rows
	ASSERT_THROW(C.bmmCaABbC(1, 1, 1, 1, B, A), MathException);
	C = Matrix(1, 2);
	// A.rows != C.rows || A.cols != C.cols
	ASSERT_THROW(C.bmmCaABbC(1, 1, 1, 1, A, B), MathException);
}
