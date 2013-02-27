#include <stdio.h>
#include <gtest/gtest.h>

using namespace std;

/*
 * MAIN - RUN ALL TESTS
 */
int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);
 #pragma warning( disable : 4244)
 	srand ( time(NULL) );

	return RUN_ALL_TESTS();
}