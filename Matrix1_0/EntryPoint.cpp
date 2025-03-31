#include "EntryPoint.h"
#include "Main.h"
#include <cstdlib> 
#include <iostream>
#include <ctime>
#include <cmath>

void Run()
{
	srand(time(0));

	// TEST
	try {

		int Mrows, Ncols;

		do {
			std::cout << "Input M rows : ";
			std::cin >> Mrows;
			std::cout << "Input N columns : ";
			std::cin >> Ncols;
		} while (Mrows <= 0 || Ncols<= 0);
		
		std::cout<< std::endl;

		std::cout << "Matrix A" << std::endl;

		Matrix A(Mrows, Ncols);
		A.FillRand();
		A.Print();

		std::cout << "SVD of matrix A: " << std::endl;

		Matrix X = A.Vmatrix(A);
		Matrix Y = A.Smatrix(A);
		Matrix Z = Z.Umatrix(A, X, Y);

		std::cout << "Matrix U : " << std::endl;
		Z.Print();
		std::cout << "Matrix S : " << std::endl;
		Y.Print();
		std::cout << "Matrix V : " << std::endl;
		X.Print();
		std::cout << "Checking" << std::endl;
		(Z * Y * (X.Transpose())).Print();
	}
	catch (const std::invalid_argument& e) {
		std::cerr << "Error : " << e.what() << std::endl;
	}

}
