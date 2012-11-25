#include <time.h>
#include "Matrix.h"
#include "CustomException.h"
#include <iostream>

using namespace std;
using namespace LibMatrix;

void createSimpleMatrixAndPrint();
void fillMatrixRandom(Matrix &m); 
void matrixMultiplication();
void checkError();
void initializeMatricesABC(Matrix &mA, Matrix &mB, Matrix &mC);

int main (int argc, char* argv[]) 
{
	srand ( time(NULL) );
	printf("\n");
	
	createSimpleMatrixAndPrint();

	system("pause");
	
	system("cls");

	matrixMultiplication();

	system("pause");

	system("cls");

	checkError();

	system("pause");

	return 0;
}

void createSimpleMatrixAndPrint()
{
	Matrix m(3,4);
	Matrix m2(3,4);

	fillMatrixRandom(m);
	fillMatrixRandom(m2);

	printf("M1: \n");
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
	double a[] = {1,2,3,4};
	Matrix testA(a);
	testA.printm("Test Matrix A with easy Constructor...\n");

	Matrix mA(2,2);
	Matrix mB(2,2);
	Matrix mC(2,2);

	initializeMatricesABC(mA,mB,mC);

	printf("Matrixmultiplikation!!! ");

	mA.printm("\nMatrix A: \n");
	mB.printm("\nMatrix B: \n");
	mC.printm("\nMatrix C: \n");

	printf("\n\nMatrix C' = alpha*A*B + beta*C: \n");
	printf("mmCaABbC\n");
	mC.mmCaABbC(3,4,mA,mB);
	mC.printm("\nMatrix C' = 3*A*B + 4*C: \n");

	initializeMatricesABC(mA,mB,mC);
	printf("\n\nMatrix C' = alpha*A' *A + beta*C: \n");
	printf("mmCaAATbC\n");
	mC.mmCaAATbC(2,3,mA);
	mC.printm("");

	initializeMatricesABC(mA,mB,mC);

	printf("\n\nMatrix C' = alpha*A' *A + beta*C: \n");
	printf("mmCaATAbC\n");
	mC.mmCaATAbC(2,3,mA);
	mC.printm("\nMatrix C' = 3*A' *A + 4*C: \n");

	initializeMatricesABC(mA,mB,mC);

	printf("\nmmCaATBbC\n");
	mC.mmCaATBbC(2,3,mA, mB);
	mC.printm("");

	initializeMatricesABC(mA,mB,mC);

	printf("\nmmCaABTbC\n");
	mC.mmCaABTbC(2,3,mA, mB);
	mC.printm("");

	initializeMatricesABC(mA,mB,mC);
	mA = mA * mB;

	mA.printm("\nMatrix A * Matrix B: \n");

	mB *= 3;
	mB.printm("\nMatrix B * 3: \n");

}

void initializeMatricesABC(Matrix &mA, Matrix &mB, Matrix &mC)
{
	double a[] = {1,2,3,4};
	double b[] = {2,5,1,3};
	double c[] = {2,3,4,1};
	mA = Matrix(a);
	mB = Matrix(b);
	mC = Matrix(c);
	// double ** data = mA.data();
	/*
	data[0][0] = 1;
	data[0][1] = 2;
	data[1][0] = 3;
	data[1][1] = 4;
	data = mB.data();
	data[0][0] = 2;
	data[0][1] = 5;
	data[1][0] = 1;
	data[1][1] = 3;
	data = mC.data();
	data[0][0] = 2;
	data[0][1] = 3;
	data[1][0] = 4;
	data[1][1] = 1;
	*/
}

void checkError() 
{
	Matrix mA(2,2);
	Matrix mB(3,3);

	try 
	{
		mA.printm("Matrix A: \n");
		mB.printm("Matrix B: \n");
		printf("\nMatrix A + Matrix B:\n");
		mA + mB;
	}
	catch (CustomException e) 
	{
		// cout << e.what();
		printf(e.what().c_str());
	}
	try 
	{
		printf("\nMatrix A - Matrix B:\n");
		mA - mB;
	}
	catch (CustomException e) 
	{
		// cout << e.what();
		printf(e.what().c_str());
	}
}

void fillMatrixRandom(Matrix &m) 
{
	double **data = m.data();
	for (int i = 0; i < m.nrows(); i++)
		for (int j = 0; j < m.ncols(); j++)
			data[i][j] = rand() % 10;
}

