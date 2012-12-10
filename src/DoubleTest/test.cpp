#include <stdio.h>
#include <gtest/gtest.h>

#include "MathException.h"
#include "Matrix.h"
using namespace LibMatrix;

using namespace std;

/*
 * MAIN - RUN ALL TESTS
 */
int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);
 
 	srand ( time(NULL) );
 
	double a[] = {
		1, 2,
		3, 4
	};
	Matrix A(2, 2, a);
 	
 	try
 	{
 		double works = A(0, 1);
 		double e = A(-1, 1);
 	}
 	catch (MathException me)
 	{
 		printf("MathException.what() = %s\n", me.what());
 		return -1;
 	}

	return RUN_ALL_TESTS();
}