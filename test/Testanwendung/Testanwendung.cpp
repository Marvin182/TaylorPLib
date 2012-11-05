#include <time.h>
#include "Matrix.h"

using namespace std;
using namespace LibMatrix;

void createSimpleMatrixAndPrint();
void fillMatrixRandom(Matrix &m); 
void matrixMultiplication();

int main (int argc, char* argv[]) 
{
	srand ( time(NULL) );
	printf("\n");
	// createSimpleMatrixAndPrint();

	// system("pause");

	system("cls");

	matrixMultiplication();

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

	printf("\n\nMatrix 1 * 2 : \n");
	m *= 2;

	m.printm("");

	printf("\n\nMatrix 3 = m * 2: \n");
	Matrix m3 = m*2;
	m3.printm("");

	printf("\n\nMatrix 4 = m4(m3): \n");
	Matrix m4(m3);
	m4.printm("");

	printf("\n\nMatrix 4 = Matrix 3 ??: \n");
	if (m4 == m3)
		printf("Are equal...");
	else
		printf("Not equal...");

	printf("\n\nMatrix 2 = Matrix 3 ??: \n");
	if (m2 == m3)
		printf("Are equal...");
	else
		printf("Not equal...");

	printf("\n\nMatrix 4 != Matrix 3 ??: \n");
	if (m4 != m3)
		printf("Not equal...");
	else
		printf("Are equal...");

	printf("\n\nMatrix 2 != Matrix 3 ??: \n");
	if (m2 != m3)
		printf("Not equal...");
	else
		printf("Are equal...");

	printf("\n\n (-) Matrix 4: \n");
	Matrix m5(-m4);
	m5.printm("");

}

void matrixMultiplication()
{
	Matrix mA(2,2);
	Matrix mB(2,2);
	double ** data = mA.data();
	data[0][0] = 1;
	data[0][1] = 2;
	data[1][0] = 3;
	data[1][1] = 4;
	data = mB.data();
	data[0][0] = 2;
	data[0][1] = 5;
	data[1][0] = 1;
	data[1][1] = 3;

	printf("Matrixmultiplikation!!! ");

	printf("\n\nMatrix A: \n");
	mA.printm("");
	printf("\n\nMatrix B: \n");
	mB.printm("");
 }

void fillMatrixRandom(Matrix &m) 
{
	double **data = m.data();
	for (int i = 0; i < m.nrows(); i++)
		for (int j = 0; j < m.ncols(); j++)
			data[i][j] = rand() % 10;
}

