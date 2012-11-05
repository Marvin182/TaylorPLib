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
}
