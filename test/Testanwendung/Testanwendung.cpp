#include <time.h>
#include "Matrix.h"

using namespace std;
using namespace LibMatrix;

void createSimpleMatrixAndPrint();
void fillMatrixRandom(Matrix &m); 

int main (int argc, char* argv[]) 
{
	srand ( time(NULL) );
	printf("\n");
	createSimpleMatrixAndPrint();

	system("pause");
	return 0;
}

void createSimpleMatrixAndPrint()
{
	Matrix m(3,4);
	Matrix m2(3,4);

	fillMatrixRandom(m);
	fillMatrixRandom(m2);

	printf("Matrix 1: \n");
	m.printm("");
	printf("\n\nMatrix 2: \n");
	m2.printm("");

	m+=m2;

	printf("\n\nMatrix 1 + Matrix 2: \n");
	m.printm("");

	printf("\n\nAddierte Matrix - Matrix 2 (= Matrix1): \n");
	m-=m2;
	m.printm("");
}

void fillMatrixRandom(Matrix &m) 
{
	double **data = m.data();
	for (int i = 0; i < m.nrows(); i++)
		for (int j = 0; j < m.ncols(); j++)
			data[i][j] = rand() % 10;
}
