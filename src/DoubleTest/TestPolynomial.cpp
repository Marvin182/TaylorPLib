#include <stdio.h>
#include <gtest/gtest.h>
#include "Polynomial.h"

using namespace std;
using namespace LibMatrix;

namespace LibMatrix {
	// needed by GTest to print matrices if some ASSERT_EQ(Polynomial1, Polynomial2) failed
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

TEST(PolynomialConstructor, default_constructor)
{
	// like Polynomial (const=1,order=0,nrcoeffs=order + 1);
	Polynomial p;
	ASSERT_EQ(true, p.isConst());
	ASSERT_EQ(0,p.order());
	ASSERT_EQ(1,p.ncoeff());
}

TEST(PolynomialConstructor, order_constructor)
{
	Polynomial p(2);
	// should also be constant, because all values are zero initialized
	ASSERT_EQ(true, p.isConst());
	ASSERT_EQ(2,p.order());
	ASSERT_EQ(3,p.ncoeff());
}

TEST(PolynomialConstructor, test_and_copy_constructor)
{
	double c[] = { 2, 3, 4 };
	Polynomial p(2);
	p.setCoeffs(c);

	ASSERT_EQ(2, p[0]);
	ASSERT_EQ(3, p[1]);
	ASSERT_EQ(4, p[2]);

	// Polynomial clone(p);
	// ASSERT_EQ(p,clone);
}
