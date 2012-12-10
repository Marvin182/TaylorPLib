#include <stdio.h>
#include <gtest/gtest.h>
#include "Polynomial.h"

using namespace std;
using namespace LibMatrix;

namespace LibMatrix {
	// needed by GTest to print polynomials if some ASSERT_EQ(Polynomial1, Polynomial2) failed
	std::ostream& operator<<(std::ostream &out, Polynomial &p)
	{
		out << setiosflags(ios::fixed) << setprecision(2);

		for (int i = 0; i < p.ncoeff(); i++)
		{
			out << p[i] << '\t';
		}
		return out; 
	}
}

/*
 * FIXTURES
 */
class PolynomialOperator: public ::testing::Test
{
protected:
	Polynomial A, B, AplusB, minusA, AminusB, twoA, AB;
	PolynomialOperator()
	{
		A = Polynomial(1);
		B = Polynomial(1);
		AplusB = Polynomial(1);
		minusA = Polynomial(1);
		AminusB = Polynomial(1);
		twoA = Polynomial(1);
		AB = Polynomial(1);

		// A
		double a[] = { 1, 2 };
		A.setCoeffs(a);

		// B
		double b[] = { 5, 2 };
		B.setCoeffs(b);

		// A + B
		double aplusb[] = { 6, 4 };
		AplusB.setCoeffs(aplusb);

		// - A
		double minusa[] = { -1, -2 };
		minusA.setCoeffs(minusa);

		// A - B
		double aminusb[] = { -4, 0 };
		AminusB.setCoeffs(aminusb);

		// 2 * A
		double twoa[] = { 2, 4 };
		twoA.setCoeffs(twoa);

		// A * B
		// echtes A * B
		// double ab[] = { 5, 12, 4 };
		// AB.setCoeffs(ab);
		// unechtes A * B (grad bleibt erhalten)
		double ab[] = { 5, 12 };
		AB.setCoeffs(ab);
	}
};

/*
 * Constructor Tests
 */
TEST(PolynomialConstructor, default_constructor)
{
	// like Polynomial (const=1,order=0,nrcoeffs=order + 1);
	Polynomial p;
	ASSERT_TRUE(p.isConst());
	ASSERT_EQ(0,p.order());
	ASSERT_EQ(1,p.ncoeff());
}

TEST(PolynomialConstructor, order_constructor)
{
	Polynomial p(2);
	// should also be constant, because all values are zero initialized
	ASSERT_TRUE(p.isConst());
	ASSERT_EQ(2,p.order());
	ASSERT_EQ(3,p.ncoeff());
}

TEST(PolynomialConstructor, test_and_copy_constructor)
{
	double c[] = { 2, 3, 4 };
	Polynomial p(2);
	p.setCoeffs(c);

	ASSERT_FALSE(p.isConst());
	ASSERT_EQ(2, p[0]);
	ASSERT_EQ(3, p[1]);
	ASSERT_EQ(4, p[2]);

	Polynomial clone(p);
	ASSERT_EQ(p,clone);
}

/*
 * Operator Tests
 */
TEST_F(PolynomialOperator, assignment)
{
	B = A;
	ASSERT_EQ(A,B);

	for (int i = 0; i < A.ncoeff(); i++)
	{
		ASSERT_EQ(A[i], B[i]);
	}
}

TEST_F(PolynomialOperator, comparsion)
{
	ASSERT_TRUE(A == A);	
	ASSERT_FALSE(A == B);

	ASSERT_FALSE(A != A);
	ASSERT_TRUE(A != B);

	Polynomial asA(A);
	ASSERT_TRUE(A == asA);

	Polynomial biggerA(2);
	for (int i = 0; i < A.ncoeff(); i++)
	{
		biggerA[i] = A[i];
	}
	// biggerA.print();
	ASSERT_FALSE(A == biggerA);

	// ASSERT_EQ and ASSERT_EQ should use the comparision operators
	ASSERT_EQ(A, asA);
	ASSERT_NE(A, B);

	ASSERT_TRUE( A < B );
	ASSERT_FALSE( A < asA );
	ASSERT_TRUE( A <= asA );
	ASSERT_TRUE( B > A );
	ASSERT_FALSE( asA > A );
	ASSERT_TRUE( asA >= A );
}

TEST_F(PolynomialOperator, plus)
{
	ASSERT_EQ(AplusB, A + B);
	A += B;
	ASSERT_EQ(AplusB, A);
}

TEST_F(PolynomialOperator, unary)
{
	ASSERT_EQ(minusA, -A);
}

TEST_F(PolynomialOperator, minus)
{
	ASSERT_EQ(AminusB, A - B);

	A -= B;
	ASSERT_EQ(AminusB, A);
}

TEST_F(PolynomialOperator, timesScalar)
{
	ASSERT_EQ(twoA, A * 2.0);
	A *= 2.0;
	ASSERT_EQ(twoA, A);
}

TEST_F(PolynomialOperator, timesPolynomial)
{
	Polynomial P1(2);
	Polynomial P2(2);
	// Polynomial PExpect(4);
	Polynomial PExpect(2);
	double x1[] = {1,2,3};
	double x2[] = {2,3,4};
	//double expect[] = {2, 7, 16, 17, 12};
	double expect[] = {2, 7, 16};
	/*
	printf("\nP1: \n");
	P1.print();
	printf("\nP2: \n");
	P2.print();
	printf("\nP1 * P2: \n");
	(P1 * P2).print();
	printf("\n");
	*/
	P1.setCoeffs(x1);
	P2.setCoeffs(x2);
	PExpect.setCoeffs(expect);
	ASSERT_EQ(PExpect, P1 * P2);

	A *= B;
	ASSERT_EQ(AB, A);
	/*

	double testX1[] = {2,3,4,5};
	double testX2[] = {4,5,6,7};

	Polynomial TestX1(3);
	Polynomial TestX2(3);
	TestX1.setCoeffs(testX1);
	TestX2.setCoeffs(testX2);
	*/
	/*
	printf("\nA: \n");
	TestX1.print();
	printf("\nB: \n");
	TestX2.print();
	printf("\nA * B: \n");
	(TestX1 * TestX2).print();
	printf("\n");
	*/
}

TEST_F(PolynomialOperator, divPolynomial)
{
	Polynomial P1(2);
	Polynomial P2(2);

	double x1[] = {2,4,6};
	double x2[] = {1,2,3};

	P1.setCoeffs(x1);
	P2.setCoeffs(x2);

	Polynomial P3 = P1 / P2;
	Polynomial P4 = P2 * P3;
	ASSERT_EQ(P1, P4);

	/*
	printf("\nP1: \n");
	P1.print();
	printf("\nP2: \n");
	P2.print();
	printf("\nP3 = P1 / P2: \n");
	P3.print();
	printf("\nP4 = P1 = P3 * P2: \n");
	P4.print();
	printf("\n");
	*/
	P1 /= P2;
	Polynomial PExpect(2);
	double expect[] = {2,0,0};
	PExpect.setCoeffs(expect);
	ASSERT_EQ(PExpect, P1);
}

/*
 * Function Tests
 */
