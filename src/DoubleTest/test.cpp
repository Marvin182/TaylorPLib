#include <stdio.h>
#include <gtest/gtest.h>
#include "Matrix.h"

using namespace std;
using namespace LibMatrix;

namespace LibMatrix {
	// needed by GTest to print matrices if some ASSERT_EQ(matrix1, matrix2) failed
	std::ostream& operator<<(std::ostream &out, const Matrix &m) {
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
		Matrix A, B, C;
		double alpha, beta;

		double getRandomNumber() { return rand() % 10; }
		void fillWithRandoms(Matrix &m) {
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
			A = Matrix(6, 6);
			fillWithRandoms(A);

			B = Matrix(6, 6);
			fillWithRandoms(B);

			C = Matrix(6, 6);
			fillWithRandoms(C);
		}
};

/*
 * CONSTRUCTOR TESTS
 */
TEST(MatrixConstrutor, default_constructor) {
	// default constructor should create a zero initialized 1x1 matrix 
	Matrix d;
	ASSERT_EQ(1, d.nrows());
	ASSERT_EQ(1, d.ncols());
	ASSERT_EQ(0, d(0, 0));
}

TEST(MatrixConstrutor, regular_constructor) {	
	// regular constructor
	Matrix r(2, 3);
	ASSERT_EQ(2, r.nrows());
	ASSERT_EQ(3, r.ncols());
	ASSERT_EQ(0, r(0, 0));
	ASSERT_EQ(0, r(1, 2));
}

TEST(MatrixConstrutor, test_and_copy_constructor) {
	// test constructor
	double values[] = {
		1, 2,	// first row
		3, 4,	// second row
		5, 6	// third row
	};
	Matrix t(3, 2, values); // rows, columns, values
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

TEST_F(MatrixOperator, assignment) {
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

TEST_F(MatrixOperator, comparsion) {
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

TEST_F(MatrixOperator, plus) {
	ASSERT_EQ(AplusB, A + B);

	A += B;
	ASSERT_EQ(AplusB, A);
}

TEST_F(MatrixOperator, unary) {
	ASSERT_EQ(minusA, -A);
}

TEST_F(MatrixOperator, minus) {
	ASSERT_EQ(AminusB, A - B);

	A -= B;
	ASSERT_EQ(AminusB, A);
}

TEST_F(MatrixOperator, timesScalar) {
	ASSERT_EQ(twoA, A * 2.0);

	A *= 2.0;
	ASSERT_EQ(twoA, A);
}

TEST_F(MatrixOperator, timesMatrix) {
	ASSERT_EQ(AB, A * B);
	ASSERT_EQ(BA, B * A);
}

/*
 * SPECIAL MATRIX MULTIPLICATIONS
 */
TEST_F(MatrixMultiplication, mmCaABbC) {
	Matrix expect = (A * B * alpha) + (C * beta);
	C.mmCaABbC(alpha, beta, A, B);

	ASSERT_EQ(expect, C);
}


TEST_F(MatrixMultiplication, bmmCaABbC) {
	// special form for B
	double b[] = {
		1, 2, 3, 0, 0, 0,
		4, 5, 6, 0, 0, 0,
		0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 1, 0,
		0, 0, 0, 0, 0, 1
	};
	B = Matrix(6, 6, b);

	(A * B).print("A * B");
	Matrix expect = (A * B * alpha) + (C * beta);
	C.bmmCaABbC(2, 3, alpha, beta, A, B);

	ASSERT_EQ(expect, C);	
}


/*
 * MAIN - RUN ALL TESTS
 */
int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
 
 	srand ( time(NULL) );
 
	return RUN_ALL_TESTS();
}