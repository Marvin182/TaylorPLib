#include <stdio.h>
#include "TaylorPLib.h"

using namespace std;
using namespace LibMatrix;

int main(int argc, char **argv)
{
	Polynomial p(3, 3.41, 4.0, 3.0, 2.0);
	
	cout << endl << p;

	return 0;
}