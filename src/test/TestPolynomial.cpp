#include <stdio.h>
#include <gtest/gtest.h>
#include "TaylorPLib.h"

using namespace std;
using namespace LibMatrix;

/*
 * FIXTURES
 */
class PolynomialOperator: public ::testing::Test
{
protected:
	Polynomial A, B, AplusB, minusA, AminusB, twoA, halfA, AB;
	
	PolynomialOperator()
	{
		A = Polynomial(1, 1.0, 2.0);
		B = Polynomial(1, 5.0, 2.0);
		AplusB = Polynomial(1, 6.0, 4.0);
		minusA = Polynomial(1, -1.0, -2.0);
		AminusB = Polynomial(1, -4.0, 0.0);
		twoA = Polynomial(1, 2.0, 4.0);
		halfA = Polynomial(1, 0.5, 1.0);
		AB = Polynomial(1, 5.0, 12.0); // we ignore the square part
	}
};

class PolynomialMethods: public ::testing::Test
{
protected:
	Polynomial ID, Zero;
	PolynomialMethods()
	{
		ID = Polynomial(0, 1.0);
		Zero = Polynomial(2, 0.0, 0.0, 0.0);
	}
};

class PolynomialExceptions: public ::testing::Test
{
protected:
	Polynomial A, B;
	PolynomialExceptions()
	{
		A = Polynomial(2, 0.0, 1.0, 2.0);
		B = Polynomial(3, 2.0, 1.0, 2.0, 4.0);
	}
};


/*
 * Constructor Tests
 */
TEST(PolynomialConstructor, default_constructor)
{
	// like Polynomial (const=1, order=0);
	Polynomial p;
	ASSERT_TRUE(p.isConst());
	ASSERT_EQ(0, p.order());
	ASSERT_EQ(1, p[0]);
}

TEST(PolynomialConstructor, order_constructor)
{
	Polynomial p(2);
	// should also be constant, because all values are zero initialized
	ASSERT_TRUE(p.isConst());
	ASSERT_EQ(2, p.order());
	ASSERT_EQ(0, p[0]);
	ASSERT_EQ(0, p[1]);
	ASSERT_EQ(0, p[2]);
}

TEST(PolynomialConstructor, short_constructor)
{
	Polynomial p(2, 2.0, 3.0, 4.0);

	ASSERT_FALSE(p.isConst());
	ASSERT_EQ(2, p[0]);
	ASSERT_EQ(3, p[1]);
	ASSERT_EQ(4, p[2]);
}

TEST(PolynomialConstructor, copy_constructor)
{
	Polynomial p(2, 2.0, 3.0, 4.0);

	ASSERT_FALSE(p.isConst());
	ASSERT_EQ(2, p[0]);
	ASSERT_EQ(3, p[1]);
	ASSERT_EQ(4, p[2]);

	Polynomial copy(p);
	ASSERT_EQ(p, copy);

	ASSERT_FALSE(copy.isConst());
	ASSERT_EQ(2, copy[0]);
	ASSERT_EQ(3, copy[1]);
	ASSERT_EQ(4, copy[2]);
}

/*
 * Operator Tests
 */
TEST_F(PolynomialOperator, assignment)
{
	B = A;
	ASSERT_EQ(A, B);
	ASSERT_EQ(B, A);

	for (int i = 0; i <= A.order(); i++)
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
	for (int i = 0; i <= A.order(); i++)
	{
		biggerA[i] = A[i];
	}
	ASSERT_FALSE(A == biggerA);

	// ASSERT_EQ and ASSERT_NE should use the comparision operators
	ASSERT_EQ(A, asA);
	ASSERT_NE(A, B);

	ASSERT_TRUE(A < B);
	ASSERT_FALSE(A < asA);
	ASSERT_TRUE(A <= asA);
	ASSERT_TRUE(B > A);
	ASSERT_FALSE(asA > A);
	ASSERT_TRUE(asA >= A);
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
	Polynomial P1(2, 1.0, 2.0, 3.0);
	Polynomial P2(2, 2.0, 3.0, 4.0);
	Polynomial PExpect(2, 2.0, 7.0, 16.0);

	ASSERT_EQ(PExpect, P1 * P2);

	A *= B;
	ASSERT_EQ(AB, A);
}

TEST_F(PolynomialOperator, divScalar)
{
	ASSERT_EQ(halfA, A / 2.0);

	A /= 2.0;
	ASSERT_EQ(halfA, A);
}

TEST_F(PolynomialOperator, divPolynomial)
{
	Polynomial P1(2, 2.0, 4.0, 6.0);
	Polynomial P2(2, 1.0, 2.0, 3.0);

	Polynomial P3 = P1 / P2;
	Polynomial P4 = P2 * P3;
	ASSERT_EQ(P1, P4);

	P1 /= P2;
	Polynomial PExpect(2, 2.0, 0.0, 0.0);
	ASSERT_EQ(PExpect, P1);
}

/*
 * Function Tests
 */
TEST_F(PolynomialMethods, isConst)
{
	Polynomial P1;
	Polynomial P2(2, 5.0, 0.0, 0.0);
	Polynomial P3(2, 5.0, 3.0, 1.0);

	ASSERT_TRUE(P1.isConst());
	ASSERT_TRUE(P2.isConst());
	ASSERT_FALSE(P3.isConst());

	Polynomial P4 = P3;
	ASSERT_FALSE(P4.isConst());
	P4 -= Polynomial(2, 0.0, 3.0, 1.0);
	ASSERT_TRUE(P4.isConst());
	P4 += Polynomial(2, 2.0, 0.0, 0.0);
	ASSERT_TRUE(P4.isConst());
	P4 -= Polynomial(2, 2.0, 1.0, 0.0);
	ASSERT_FALSE(P4.isConst());
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
	P1.set2Zero();
	ASSERT_TRUE(P1.isZero(1));
	P1.setCoeffs(p1);
	ASSERT_FALSE(P1.isZero(1));
	P1.set2Zero(2);
	ASSERT_TRUE(P1.isZero(1));
}

TEST_F(PolynomialMethods, set2Const)
{
	Polynomial P1(3);
	double p1[] = {1,1,2,1};
	P1.setCoeffs(p1);

	ASSERT_FALSE(P1.isConst());
	P1.set2Const(6);
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