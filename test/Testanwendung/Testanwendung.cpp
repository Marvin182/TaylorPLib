#include "Matrix.h"

using namespace std;
using namespace LibMatrix;

void createSimpleMatrixAndPrint();

int main (int argc, char* argv[]) 
{
	printf("\n");
	createSimpleMatrixAndPrint();

	system("pause");
	return 0;
}

void createSimpleMatrixAndPrint()
{
	Matrix m(3,4);
	m.printm("");
	double ** data = m.data();
	data[0][0] = 15;

	m(0,1) = 12;

	m.printm("");
	printf("m[0][0] = %f\n", data[0][0]);
}
