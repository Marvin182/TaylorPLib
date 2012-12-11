#include <stdio.h>
#include <gtest/gtest.h>
#include "TaylorPLib.h"

using namespace std;
using namespace LibMatrix;

namespace LibMatrix {
	// needed by GTest to print polynomials if some ASSERT_EQ(Polynomial1, Polynomial2) failed
	std::ostream& operator<<(std::ostream &out, const Polynomial &p)
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

class PolynomialMethods: public ::testing::Test
{
protected:
	Polynomial ID, Zero;
	PolynomialMethods() {

		Zero = Polynomial(2);

		// ID 
		double id[] = {1};
		ID.setCoeffs(id);

		// Zero
		double zero[] = {0,0,0};
		Zero.setCoeffs(zero);
	}
};

class PolynomialExceptions: public ::testing::Test
{
protected:
	Polynomial A, B;
	PolynomialExceptions()
	{
		double a[] = {0,1,2};
		double b[] = {2,1,2,4};
		A = Polynomial(2);
		B = Polynomial(3);
		A.setCoeffs(a);
		B.setCoeffs(b);
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
TEST_F(PolynomialMethods, isConst)
{
	Polynomial P1;
	Polynomial P2(2);
	Polynomial P3(2);

	double p2[] = {5,0,0};
	P2.setCoeffs(p2);

	double p3[] = {5,3,1};
	P3.setCoeffs(p3);

	ASSERT_TRUE(P1.isConst());
	ASSERT_TRUE(P2.isConst());
	ASSERT_FALSE(P3.isConst());
}

TEST_F(PolynomialMethods, isConstWithEps)
{
	Polynomial P1(2);
	Polynomial P2(2);

	double p1[] = {4,2,3};
	P1.setCoeffs(p1);

	double p2[] = {5,3,1};
	P2.setCoeffs(p2);

	ASSERT_TRUE(P1.isConst(3));
	ASSERT_FALSE(P2.isConst(2));
}

TEST_F(PolynomialMethods, isId)
{
	ASSERT_TRUE(ID.isId());
	Polynomial P1(4);
	Polynomial P2(3);
	double p1[] = {1,0,0,0,0};
	P1.setCoeffs(p1);
	double p2[] = {1,2,1,3,2};
	P2.setCoeffs(p2);

	ASSERT_TRUE(P1.isId());
	ASSERT_FALSE(P2.isId());
}

TEST_F(PolynomialMethods, isIdWithEps)
{
	ASSERT_TRUE(ID.isId());
	Polynomial P1(4);
	double p1[] = {1,2,3,4,3};
	P1.setCoeffs(p1);

	ASSERT_TRUE(P1.isId(4));
	ASSERT_FALSE(P1.isId(3));
}

TEST_F(PolynomialMethods, sqr)
{
	Polynomial P1(3);
	double p1[] = {1,2,3,4};
	P1.setCoeffs(p1);
	Polynomial P2 = P1.sqr();
	Polynomial PExpect(3);
	double pexpect[] = {1,4,10,20};
	PExpect.setCoeffs(pexpect);
	ASSERT_EQ(P2,PExpect);
	double pnotexpect[] = {1,4,10,30};
	PExpect.setCoeffs(pnotexpect);
	ASSERT_NE(P2,PExpect);
}

TEST_F(PolynomialMethods, setsqr)
{
	Polynomial P1(3);
	double p1[] = {1,2,3,4};
	P1.setCoeffs(p1);
	Polynomial P2 = P1.sqr();
	P1.setSqr();
	ASSERT_EQ(P1,P2);
}

TEST_F(PolynomialMethods, sqrt)
{
	Polynomial P1(3);
	double p1[] = {1,4,10,20};
	P1.setCoeffs(p1);
	Polynomial P2 = P1.sqrt();
	Polynomial PExpect(3);
	double pexpect[] = {1,2,3,4};
	PExpect.setCoeffs(pexpect);
	ASSERT_EQ(P2,PExpect);

	double pnotexpect[] = {1,2,3,5};
	PExpect.setCoeffs(pnotexpect);
	ASSERT_NE(P2,PExpect);
}

TEST_F(PolynomialMethods, setsqrt)
{
	Polynomial P1(3);
	double p1[] = {1,4,10,20};
	P1.setCoeffs(p1);
	Polynomial P2 = P1.sqrt();
	P1.setSqrt();
	ASSERT_EQ(P1,P2);
}

TEST_F(PolynomialMethods, isZero)
{
	ASSERT_TRUE(Zero.isZero());
}

TEST_F(PolynomialMethods, isZeroWithEps)
{
	Polynomial P1(3);
	double p1[] = {1,1,2,1};
	P1.setCoeffs(p1);

	ASSERT_TRUE(P1.isZero(2));
	ASSERT_FALSE(P1.isZero(1));
}

TEST_F(PolynomialMethods, set2Zero)
{
	Polynomial P1(3);
	double p1[] = {1,1,2,1};
	P1.setCoeffs(p1);
	ASSERT_FALSE(P1.isZero(1));
	P1.set2zero();
	ASSERT_TRUE(P1.isZero(1));
	P1.setCoeffs(p1);
	ASSERT_FALSE(P1.isZero(1));
	P1.set2zero(2);
	ASSERT_TRUE(P1.isZero(1));
}

TEST_F(PolynomialMethods, set2Const)
{
	Polynomial P1(3);
	double p1[] = {1,1,2,1};
	P1.setCoeffs(p1);

	ASSERT_FALSE(P1.isConst());
	P1.set2const(6);
	ASSERT_TRUE(P1.isConst());

	Polynomial P2(3);
	double p2[] = {6,0,0,0};
	P2.setCoeffs(p2);

	ASSERT_EQ(P1, P2);
}

TEST_F(PolynomialMethods, feval)
{
	Polynomial P1(3);
	double p1[] = {4,3,2,1};
	P1.setCoeffs(p1);

	ASSERT_EQ(4,P1.feval());
}

TEST_F(PolynomialMethods, eval)
{
	Polynomial P1(3);
	double p1[] = {1,2,3,4};
	P1.setCoeffs(p1);

	ASSERT_EQ(10,P1.eval(2,1));
	ASSERT_NE(10,P1.eval(2,0.5));
}

TEST_F(PolynomialMethods, shift)
{
	Polynomial P1(3);
	double p1[] = {1,2,3,4};
	P1.setCoeffs(p1);

	Polynomial P2(3);
	double p2[] = {2,6,12,0};
	P2.setCoeffs(p2);

	P1.shift();
	ASSERT_EQ(P1, P2);

	Polynomial P3(3);
	double p3[] = {6,24,0,0};
	P3.setCoeffs(p3);
	P1.shift();
	ASSERT_EQ(P1, P3);

}

TEST_F(PolynomialExceptions, OperatorExceptions)
{
	ASSERT_THROW(A[-1],MathException);
	ASSERT_THROW(A[5],MathException);
	ASSERT_THROW(A + B,MathException);
	ASSERT_THROW(A += B,MathException);
	ASSERT_THROW(A - B,MathException);
	ASSERT_THROW(A -= B,MathException);
	ASSERT_THROW(A * B,MathException);
	ASSERT_THROW(A *= B,MathException);
	ASSERT_THROW(A / B,MathException);
	ASSERT_THROW(A /= B,MathException);
}

TEST_F(PolynomialExceptions, FunctionExceptions)
{
	ASSERT_THROW(Polynomial(-1),MathException);
}