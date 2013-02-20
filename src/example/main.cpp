#include <stdio.h>
#include "TaylorPLib.h"

using namespace std;
using namespace LibMatrix;

void showPolynom();
void showMatrix();
void showSpecialMatrixMultiplication();
void showError();

int main(int argc, char **argv)
{
	cout << "Polynome: " << endl;
	showPolynom();
	system("pause");
	system("cls");
	cout << "Matrizen mit Polynomen: " << endl;
	showMatrix();
	system("pause");
	system("cls");
	cout << "Spezielle Matrizenmultiplikation: " << endl;
	showSpecialMatrixMultiplication();
	system("pause");
	system("cls");
	cout << "Exception Catching: " << endl;
	showError();
	system("pause");

	return 0;
}

void showPolynom()
{
	Polynomial P1(2, 2.0, 4.0, 3.0);
	Polynomial P2(2, 2.0, 4.0, 3.0);

	Polynomial P3 = P1 + P2;
	Polynomial P4 = P1 * P2;
	Polynomial P5 = P4 / P2;

	cout << "P1: " << endl << P1 << endl;
	cout << "P2: " << endl << P2 << endl;
	cout << "P1 + P2: " << endl << P3 << endl;
	cout << "P1 * P2: " << endl << P4 << endl;
	cout << "P1 * P2 / P2: " << endl << P5 << endl;
}

void showMatrix()
{
	Polynomial Pa[] = {
		Polynomial(1, 1.0, 2.0), Polynomial(1, 3.0, 4.0), Polynomial(1, 5.0, 6.0), Polynomial(1, 1.0, 2.0), Polynomial(1, 5.0, 6.0),
		Polynomial(1, 1.0, 2.0), Polynomial(1, 3.0, 4.0), Polynomial(1, 5.0, 6.0), Polynomial(1, 1.0, 2.0), Polynomial(1, 5.0, 6.0),
		Polynomial(1, 1.0, 2.0), Polynomial(1, 3.0, 4.0), Polynomial(1, 5.0, 6.0), Polynomial(1, 1.0, 2.0), Polynomial(1, 5.0, 6.0),
		Polynomial(1, 1.0, 2.0), Polynomial(1, 3.0, 4.0), Polynomial(1, 5.0, 6.0), Polynomial(1, 1.0, 2.0), Polynomial(1, 5.0, 6.0),
		Polynomial(1, 1.0, 2.0), Polynomial(1, 3.0, 4.0), Polynomial(1, 5.0, 6.0), Polynomial(1, 1.0, 2.0), Polynomial(1, 5.0, 6.0)
	};
	Polynomial Pb[] = {
		Polynomial(1, 0.0, 0.0), Polynomial(1, 0.0, 0.0), Polynomial(1, 0.0, 0.0), Polynomial(1, 2.0, 1.0), Polynomial(1, 2.0, 1.0), 
		Polynomial(1, 0.0, 0.0), Polynomial(1, 0.0, 0.0), Polynomial(1, 0.0, 0.0), Polynomial(1, 2.0, 1.0), Polynomial(1, 2.0, 1.0), 
		Polynomial(1, 0.0, 0.0), Polynomial(1, 0.0, 0.0), Polynomial(1, 0.0, 0.0), Polynomial(1, 2.0, 1.0), Polynomial(1, 2.0, 1.0), 
		Polynomial(1, 0.0, 0.0), Polynomial(1, 0.0, 0.0), Polynomial(1, 0.0, 0.0), Polynomial(1, 2.0, 1.0), Polynomial(1, 2.0, 1.0), 
		Polynomial(1, 0.0, 0.0), Polynomial(1, 0.0, 0.0), Polynomial(1, 0.0, 0.0), Polynomial(1, 2.0, 1.0), Polynomial(1, 2.0, 1.0)
	};

	Matrix M1 = Matrix(5,5,Pa);
	Matrix M2 = Matrix(5,5,Pb);

	Matrix M3 = M1 + M2;
	Matrix M4 = M1 * M2;

	cout << "M1: " << endl << M1 << endl;
	cout << "M2: " << endl << M2 << endl;
	cout << "M1 + M2: " << endl << M3 << endl;
	cout << "M1 * M2: " << endl << M4 << endl;
}

void showSpecialMatrixMultiplication()
{
	Polynomial Pa[] = {
		Polynomial(1, 1.0, 2.0), Polynomial(1, 3.0, 4.0), Polynomial(1, 5.0, 6.0),
		Polynomial(1, 7.0, 8.0), Polynomial(1, 9.0, 0.0), Polynomial(1, 1.0, 2.0),
		Polynomial(1, 3.0, 4.0), Polynomial(1, 5.0, 6.0), Polynomial(1, 7.0, 8.0)
	};
	Polynomial Pb[] = {
		Polynomial(1, 3.0, 4.0), Polynomial(1, 5.0, 6.0), Polynomial(1, 7.0, 8.0),
		Polynomial(1, 1.0, 2.0), Polynomial(1, 3.0, 4.0), Polynomial(1, 5.0, 6.0),
		Polynomial(1, 7.0, 8.0), Polynomial(1, 9.0, 0.0), Polynomial(1, 1.0, 2.0)
	};
	Polynomial Pc[] = {
		Polynomial(1, 7.0, 8.0), Polynomial(1, 9.0, 0.0), Polynomial(1, 1.0, 2.0),
		Polynomial(1, 3.0, 4.0), Polynomial(1, 5.0, 6.0), Polynomial(1, 7.0, 8.0),
		Polynomial(1, 1.0, 2.0), Polynomial(1, 3.0, 4.0), Polynomial(1, 5.0, 6.0)
	};

	Matrix  M1 = Matrix(3,3,Pa);
	Matrix  M2 = Matrix(3,3,Pb);
	Matrix  M3 = Matrix(3,3,Pc);

	cout << "( M1 * M2 * alpha ) + ( M3 * beta ) where: " << endl;
	cout << "M1: " << endl << M1 << endl;
	cout << "M2: " << endl << M2 << endl;
	cout << "M3: " << endl << M3 << endl;
	cout << "alpha: 2" << endl;
	cout << "beta: 2" << endl;

	Matrix calculated = (M1 * M2 * 2) + ( M3 * 2);

	cout << "Result Normal Mode: " << endl << calculated << endl;
	M3.mmCaABbC(2.0,2.0,M1, M2);
	cout << "Result Function Mode: " << endl << M3 << endl;

}

void showError()
{
	Polynomial Pa[] = {
		Polynomial(1, 1.0, 2.0), Polynomial(1, 3.0, 4.0), Polynomial(1, 5.0, 6.0),
		Polynomial(1, 7.0, 8.0), Polynomial(1, 9.0, 0.0), Polynomial(1, 1.0, 2.0),
		Polynomial(1, 3.0, 4.0), Polynomial(1, 5.0, 6.0), Polynomial(1, 7.0, 8.0)
	};
	Polynomial Pb[] = {
		Polynomial(1, 3.0, 4.0), Polynomial(1, 5.0, 6.0),
		Polynomial(1, 1.0, 2.0), Polynomial(1, 3.0, 4.0)
	};

	Matrix M1 = Matrix(3,3,Pa);
	Matrix M2 = Matrix(2,2,Pb);

	cout << "M1: " << endl << M1 << endl;
	cout << "M2: " << endl << M2 << endl;
	cout << "M1 + M2: " << endl;

	try 
	{
		Matrix M3 = M1 + M2;
		cout << "M3: " << endl << M3 << endl; 
	} 
	catch (MathException me)
	{
		cout << "Es ist ein Fehler aufgetreten: " << endl << me.what() << endl;
	}

}
