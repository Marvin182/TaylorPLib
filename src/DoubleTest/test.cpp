#include <stdio.h>
#include <gtest/gtest.h>
#include "Matrix.h"

using namespace std;
using namespace LibMatrix;

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
			double a[] = {
				1, 2,
				3, 4
			};
			A = Matrix(2, 2, a);

			// B
			double b[] = {
				5, 2,
				0, 3
			};
			B = Matrix(2, 2, b);

			// A + B
			double aplusb[] = {
				6, 4,
				3, 7 
			};
			AplusB = Matrix(2, 2, aplusb);

			// - A
			double minusa[] = {
				-1, -2,
				-3, -4
			};
			minusA = Matrix(2, 2, minusa);

			// A - B
			double aminusb[] = {
				-4, 0,
				 3, 1
			};
			AminusB = Matrix(2, 2, aminusb);

			// 2 * A
			double twoa[] = {
				2, 4,
				6, 8
			};
			twoA = Matrix(2, 2, twoa);

			// A * B
			double ab[] = {
				 5,  8,
				15, 18
			};
			AB = Matrix(2, 2, ab);			
			
			// B * A
			double ba[] = {
				11, 18,
				 9, 12
			};
			BA = Matrix(2, 2, ba);
		}
};

class MatrixMultiplication: public ::testing::Test
{
	protected:
		Matrix A, B, C, I;
		double alpha, beta;

		double getRandomNumber() { return rand() % 10; }
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
			double i[] = {
				1, 0, 0,
				0, 1, 0,
				0, 0, 1
			};
			I = Matrix(3, 3, i);

			double z[] = {
				0, 0, 0,
				0, 0, 0,
				0, 0, 0
			};
			Z = Matrix(3, 3, z);

			double a[] = {
				1, 2, 3,
				4, 5, 6,
				7, 8, 9
			};
			A = Matrix(3, 3, a);

			double b[] = {
				1,  2,  3,  4,
				5,  6,  7,  8,
				9, 10, 11, 12
			};
			B = Matrix(3, 4, b);

			double c[] = {
				5, 5, 5, 5, 5, 5, 5,
				5, 5, 5, 5, 5, 5, 5,
				5, 5, 5, 5, 5, 5, 5,
				5, 5, 5, 5, 5, 5, 5,
				5, 5, 5, 5, 5, 5, 5,
				5, 5, 5, 5, 5, 5, 5,
				5, 5, 5, 5, 5, 5, 5,
				5, 5, 5, 5, 5, 5, 5
			};
			C = Matrix(8, 7, c);

			double cz[] = {
				5, 5, 5, 5, 5, 5, 5,
				5, 5, 0, 0, 5, 5, 5,
				5, 5, 0, 0, 5, 5, 5,
				5, 5, 0, 0, 5, 5, 5,
				5, 5, 0, 0, 5, 5, 5,
				5, 5, 5, 5, 5, 5, 5,
				5, 5, 5, 5, 5, 5, 5,
				5, 5, 5, 5, 5, 5, 5
			};
			CZ = Matrix(8, 7, cz);

			double cid[] = {
				5, 5, 5, 5, 5, 5, 5,
				5, 5, 1, 0, 0, 0, 5,
				5, 5, 0, 1, 0, 0, 5,
				5, 5, 0, 0, 1, 0, 5,
				5, 5, 0, 0, 0, 1, 5,
				5, 5, 5, 5, 5, 5, 5,
				5, 5, 5, 5, 5, 5, 5,
				5, 5, 5, 5, 5, 5, 5
			};
			CId = Matrix(8, 7, cid);

			double c3[] = {
				5, 5, 5, 5, 5, 5, 5,
				5, 5, 5, 5, 5, 5, 5,
				5, 3, 3, 3, 3, 5, 5,
				5, 3, 3, 3, 3, 5, 5,
				5, 3, 3, 3, 3, 5, 5,
				5, 5, 5, 5, 5, 5, 5,
				5, 5, 5, 5, 5, 5, 5,
				5, 5, 5, 5, 5, 5, 5
			};
			C3 = Matrix(8, 7, c3);
		}
};

/*
 * CONSTRUCTOR TESTS
 */
TEST(MatrixConstructor, default_constructor)
{
	// default constructor should create a zero initialized 1x1 matrix 
	Matrix d;
	ASSERT_EQ(1, d.nrows());
	ASSERT_EQ(1, d.ncols());
	ASSERT_EQ(0, d(0, 0));
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
	double values[] = {
		1, 2,	// first row
		3, 4,	// second row
		5, 6	// third row
	};
	Matrix t(3, 2, values); // rows, columns, values
	t.print("t");
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
	t(0, 1) = 8;
	t(2, 0) = 9;
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
	double b[] = {
		1, 2, 3, 0, 0,
		4, 5, 6, 0, 0,
		7, 8, 9, 0, 0,
		0, 0, 0, 1, 0,
		0, 0, 0, 0, 1
	};

	B = Matrix(5, 5, b);

	Matrix expect = (A * B * alpha) + (C * beta);
	C.bmmCaABbC(3, 3, alpha, beta, A, B);

	ASSERT_EQ(expect, C);	
}

TEST_F(MatrixMultiplication, mmCasABbC)
{
   	// special form for A
	double a[] = {
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		1, 2, 3, 4, 5,
		5, 4, 3, 2, 1
	};
	A = Matrix(5, 5, a);

	Matrix expect = (A * B * alpha) + (C * beta);
	C.mmCasABbC(2, alpha, beta, A, B);
	ASSERT_EQ(expect, C);	
}

TEST_F(MatrixMultiplication, mmCaAsBbC)
{
   	// special form for B
	double b[] = {
		0, 0, 0, 1, 5,
		0, 0, 0, 2, 4,
		0, 0, 0, 3, 3,
		0, 0, 0, 4, 2,
		0, 0, 0, 5, 1
	};
	B = Matrix(5, 5, b);

	Matrix expect = (A * B * alpha) + (C * beta);
	C.mmCaAsBbC(2, alpha, beta, A, B);
	ASSERT_EQ(expect, C);	
}

TEST_F(MatrixMultiplication, mmCaAUTBPbC)
{
	double bUT[] = {
		1, 2, 3, 4, 5,
		0, 6, 7, 8, 9,
		0, 0, 1, 2, 3,
		0, 0, 0, 4, 5,
		0, 0, 0, 0, 6
	};
	Matrix B_UT = Matrix(5, 5, bUT);
	
	// special form for B
	double b[] = {
		1, 3, 2, 4, 5,
		0, 7, 6, 8, 9,
		0, 1, 0, 2, 3,
		0, 0, 0, 4, 5,
		0, 0, 0, 0, 6
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
	double a1[] = {
		1, 2, 3, 4, 5,
		6, 7, 8, 9, 1,
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0
	};
	
	double b1[] = {
		1, 2, 3, 4, 5,
		6, 7, 8, 9, 1,
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0
	};
	
	double b2[] = {
		1, 2, 3, 5, 4,
		6, 7, 8, 1, 9,
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0, 0, 0, 0
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
	double b[] = {
		1, 6, 0, 0, 0,
		2, 7, 0, 0, 0,
		3, 8, 0, 0, 0,
		4, 9, 0, 0, 0,
		5, 1, 0, 0, 0
	};

	B = Matrix(5,5, b);
	Matrix bt = B.asTranspose();
    Matrix expect = (A * bt * alpha) + (C * beta);

	C.mmCaABTbC(2, true, alpha, beta, A, B);
	ASSERT_EQ(expect, C);
}

TEST_F(MatrixMultiplication, mmCaABTbC_down)
{
	double b[] = {
		0, 0, 0, 1, 6,
		0, 0, 0, 2, 7,
		0, 0, 0, 3, 8,
		0, 0, 0, 4, 9,
		0, 0, 0, 5, 1
	};

	B = Matrix(5,5, b);
	Matrix bt = B.asTranspose();
    Matrix expect = (A * bt * alpha) + (C * beta);

	C.mmCaABTbC(2, false, alpha, beta, A, B);
	ASSERT_EQ(expect, C);
}

TEST_F(MatrixMultiplication, bmmCaABTbC)
{
	// case 1
	double a[] = {
		1, 2, 3, 0, 0,
		4, 5, 6, 0, 0,
		7, 8, 9, 0, 0,
		0, 0, 0, 1, 0,
		0, 0, 0, 0, 1
	};
	A = Matrix(5, 5, a);
	
	Matrix bt = B.asTranspose();
	Matrix expect = (A * bt * alpha) + (C * beta);
	Matrix tempC = C;

	C.bmmCaABTbC(3, 3, alpha, beta, A, B);
	ASSERT_EQ(expect, C);

	tempC.bmmCaABTbC(2, 2, alpha, beta, A, B);
	ASSERT_NE(expect, tempC);

	// case 2
	double a2[] = {
		1,2,3,0,
		4,5,6,0,
		0,0,0,1
	};
	double b2[] = {
		1,2,3,
		4,5,6,
		0,0,0,
		0,0,0
	};

	A = Matrix(3, 4, a2);
	B = Matrix(4, 3, b2);
	C = Matrix(3, 3);
	// Exception, da B transposed wird und da die Mulitplikation (3 x 4) * (3 x 4) nicht möglich ist
	ASSERT_THROW(C.bmmCaABTbC(2, 3, alpha, beta, A, B), CustomException);
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
	double b[] = {
		 4,  1,  3,  2,
		 8,  5,  7,  6,
		12,  9,  11, 10
	};
	Matrix newB = Matrix(3, 4, b);

	B.cpermutem((int*) newOrder);
	ASSERT_EQ(newB, B);	
}

TEST_F(MatrixMethods, cpermutem_trans)
{
	int newOrder[] = {1, 3, 2, 0};
	// Beim Transpose wird aus dem Pivot Vektor -> { 3, 0, 2, 1 }

	double b[] = {
		 4,  1,  3,  2,
		 8,  5,  7,  6,
		12,  9,  11, 10
	};
	Matrix newB = Matrix(3, 4, b);

	B.cpermutem((int*) newOrder, true);
	ASSERT_EQ(newB, B);	
}

TEST_F(MatrixMethods, rpermutem)
{
	int newOrder[] = {1, 2, 0};
	double b[] = {
		5,  6,  7,  8,
		9, 10, 11, 12,
		1,  2,  3,  4
	};
	Matrix newB = Matrix(3, 4, b);

	B.rpermutem((int*) newOrder);
	ASSERT_EQ(newB, B);
}

TEST_F(MatrixMethods, transpose)
{
	// test in place transpose
	double atrans[] = {
		1, 4, 7,
		2, 5, 8,
		3, 6, 9
	};
	Matrix Atranspose(3, 3, atrans);
	A.transpose();
	ASSERT_EQ(Atranspose, A);

	// test normal transpose
	double btrans[] = {
		1,  5,  9,
		2,  6, 10,
		3,  7, 11,
		4,  8, 12
	};
	Matrix Btranspose(4, 3, btrans);
	B.transpose();
	ASSERT_EQ(Btranspose, B);	
}

TEST_F(MatrixMethods, asTranspose)
{
	double atrans[] = {
		1, 4, 7,
		2, 5, 8,
		3, 6, 9
	};
	Matrix Atranspose(3, 3, atrans);
	ASSERT_EQ(Atranspose, A.asTranspose());

	double btrans[] = {
		1,  5,  9,
		2,  6, 10,
		3,  7, 11,
		4,  8, 12
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

	double zeros[] = {
		0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0
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
	double b42[] = {
		42, 42, 42, 42,
		42, 42, 42, 42,
		42, 42, 42, 42,
	};
	Matrix B42(3, 4, b42);

	B.set2Val(42.0);
	ASSERT_EQ(B42, B);
}

TEST_F(MatrixMethods, set2Val2)
{
	C.set2Val(2, 3, 1, 2, 3.0);
	ASSERT_EQ(C3, C);
}

TEST_F(MatrixMethods, set2ValFromIndices)
{
	C.set2ValFromIndices(2, 4, 1, 4, 3.0);
	ASSERT_EQ(C3, C);
}

// int Matrix::trinvm(Matrix &Am);
// int Matrix::trinvm(int r, Matrix &Am);
// bool Matrix::mcompare(Matrix &B);
// bool Matrix::mcompare(Matrix &B, double eps);
// bool Matrix::mcompare(Matrix &B, int r, double eps);

/*
 * MAIN - RUN ALL TESTS
 */
int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);
 
 	srand ( time(NULL) );
 
	return RUN_ALL_TESTS();
}