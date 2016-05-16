#include <iostream>
#include "stdafx.h"
#include "Multigrid.h"
#include <fstream>
#include <time.h>

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	// Set initial values //
	int n = INITIAL_N;
	int size = n*n;
	int max_D = MAX_D;
	int iteration = ITERATION;
	int cycle = CYCLE;
	int non_zeros = max_D*size;
	clock_t start, end;
	start = clock();

	// Initialize CSR matrix A //
	CSR A_h = CSR(non_zeros,size,size);
	build_A(&A_h,n);
	
	// initialize f, solution, error //
	double * f_h = new double[size];
	double * solution = new double[size];
	double * error = new double[size];
	build_solution(solution,n);
	build_f(f_h,n);
	
	// Print information before starting V-cycle //
	cout << " ** Algebraic Multigrid Method **" << endl<<endl;
	cout << "The initial matrix size is " << n << " by " << n << endl;
	cout << "Implement ("<<iteration<<",0) V-cycle "<<cycle<<" times"<< endl<<endl;
	cout<<endl;
	
	// Implement the first V-cycle and set up operators //
	cout << "<Implement the first V-cycle and set up operators>" <<endl;
	cout<<endl;
	cout << "A_h : Number of rows = " << size << ", Number of nonzeros = " << A_h.cnt<<endl;

	// 1st : h->2h //
	int new_size;
	CSR I_h, T_h;
	new_size = coarse(&A_h, size, &I_h, &T_h);	
		
	double * f_2h = new double[new_size];
	double * u_h = new double[size];
	for(int i=0; i<size; i++)
		u_h[i]=0;
	restrict(&A_h, f_h, size, u_h, f_2h, &T_h);
	
	non_zeros = max_D*new_size;
	CSR A_2h=CSR(non_zeros);
	coarse_A(&T_h, &A_h, &I_h, &A_2h); //  A_2h <- T_h * A_h * I_h
	cout << "A_2h : Number of rows = " << new_size << ", Number of nonzeros = " << A_2h.cnt<<endl;

	// 2nd : 2h->4h //
	size = new_size;
	CSR I_2h, T_2h;
	new_size = coarse(&A_2h, size, &I_2h, &T_2h);	
		
	double * f_4h = new double[new_size];
	double * u_2h = new double[size];
	for(int i=0; i<size; i++)
		u_2h[i]=0;
	restrict(&A_2h, f_2h, size, u_2h, f_4h, &T_2h);
	
	non_zeros = max_D*new_size;
	CSR A_4h=CSR(non_zeros);
	coarse_A(&T_2h, &A_2h, &I_2h, &A_4h); //  A_4h <- T_2h * A_2h * I_2h
	cout << "A_4h : Number of rows = " << new_size << ", Number of nonzeros = " << A_4h.cnt<<endl;
		
	// 3rd : 4h->8h //
	size = new_size;
	CSR I_4h, T_4h;
	new_size = coarse(&A_4h, size, &I_4h, &T_4h);	
		
	double * f_8h = new double[new_size];
	double * u_4h = new double[size];
	for(int i=0; i<size; i++)
		u_4h[i]=0;
	restrict(&A_4h, f_4h, size, u_4h, f_8h, &T_4h);
	
	non_zeros = max_D*new_size;
	CSR A_8h=CSR(non_zeros);
	coarse_A(&T_4h, &A_4h, &I_4h, &A_8h); //  A_8h <- T_4h * A_4h * I_4h
	cout << "A_8h : Number of rows = " << new_size << ", Number of nonzeros = " << A_8h.cnt<<endl;
	
	// 4th : 8h->16h //
	size = new_size;
	CSR I_8h, T_8h;
	new_size = coarse(&A_8h, size, &I_8h, &T_8h);	
		
	double * f_16h = new double[new_size];
	double * u_8h = new double[size];
	for(int i=0; i<size; i++)
		u_8h[i]=0;
	restrict(&A_8h, f_8h, size, u_8h, f_16h, &T_8h);
	
	non_zeros = max_D*new_size;
	CSR A_16h=CSR(non_zeros);
	coarse_A(&T_8h, &A_8h, &I_8h, &A_16h); //  A_16h <- T_8h * A_8h * I_8h
	cout << "A_16h : Number of rows = " << new_size << ", Number of nonzeros = " << A_16h.cnt<<endl;
	cout<< endl;

	// Calculate error //
	bool tracking = TRACKING;
	if(tracking)
		cout << "Calculating error...";
	size = A_16h.r;
	double * e_16h = new double[size];
	for(int i=0; i<size; i++)
		e_16h[i]=0;
	
	// Implement Gauss-Seidel with criterion on the bottom level //
	double criterion = 1/1000000;
	int iter = A_16h.gauss_seidel_criterion(e_16h,f_16h,criterion);  // Solve, A_16h * e_16h = f_16h , so we get e_16h
	if(tracking)
		cout << "Done"<<endl<<" - Implemented Gauss-Seidel " << iter <<" times with criterion 10^-6"<<endl<<endl;

	// Interpolate the coarse-grid error to the fine grid //

	// 16h->8h //
	size = A_8h.r;
	double * e_8h = new double[size];
	I_8h.matrix_by_vector(e_16h, e_8h, A_16h.r);// e_8h <- I_8h * e_16h

	// 8h->4h //
	for(int i=0; i<size; i++)
		u_8h[i] += e_8h[i];						// u_8h <- u_8h + e_8h
	size = A_4h.r;
	double * e_4h = new double[size];
	I_4h.matrix_by_vector(u_8h, e_4h, A_8h.r);	// e_4h <- I_4h * e_8h

	// 4h->2h //
	for(int i=0; i<size; i++)
		u_4h[i] += e_4h[i];						// u_4h <- u_4h + e_4h
	size = A_2h.r;
	double * e_2h = new double[size];
	I_2h.matrix_by_vector(u_4h, e_2h, A_4h.r);	 // e_2h <- I_2h * e_4h

	// 2h->h //
	for(int i=0; i<size; i++)
		u_2h[i] += e_2h[i];						// u_2h <- u_2h + e_2h
	size = A_h.r;
	double * e_h = new double[size];
	I_h.matrix_by_vector(u_2h, e_h, A_2h.r);	 // e_h <- I_h * u_2h
	
	// get u_h //
	for(int i=0; i<size; i++)
		u_h[i] += e_h[i];						// u_h <- u_h + e_h

	// Get error vector and print it //
	for(int i=0; i<size; i++)
		error[i]=solution[i]-u_h[i];
	
	cout << endl << "<Implement V-cycle "<< cycle <<" times>"<<endl;
	cout <<endl;
	cout << " 1-th V-cycle";
	cout << " : Norm of error vector = "<<norm(error,size)<<endl;
	
	// Implement V-cycle the remaining CYCLE-times //
	for(int i=0; i<cycle-1; i++)
	{
		cout << " "<< i+2 <<"-th V-cycle";
		size = A_h.r;
		build_f(f_h,n);
		A_h.gauss_seidel(u_h,f_h,iteration);
		make_f_remainder(&A_h,u_h,f_h,size); // f <- (f-Au)
		T_h.matrix_by_vector(f_h, f_2h, size); // f_2h <- I_t*f
	
		size = A_2h.r;
		for(int i=0; i<size; i++)
			u_2h[i]=0;
		A_2h.gauss_seidel(u_2h,f_2h,iteration);
		make_f_remainder(&A_2h,u_2h,f_2h,size); // f_2h <- (f_2h - A_2h * u_2h)
		T_2h.matrix_by_vector(f_2h, f_4h, size); // f_4h <- I_t_2h * f_2h
	
		size = A_4h.r;
		for(int i=0; i<size; i++)
			u_4h[i]=0;
		A_4h.gauss_seidel(u_4h,f_4h,iteration);
		make_f_remainder(&A_4h,u_4h,f_4h,size); // f_4h <- (f_4h - A_4h * u_4h)
		T_4h.matrix_by_vector(f_4h, f_8h, size); // f_8h <- I_t_4h * f_4h
	
		size = A_8h.r;
		for(int i=0; i<size; i++)
			u_8h[i]=0;
		A_8h.gauss_seidel(u_8h,f_8h,iteration);
		make_f_remainder(&A_8h,u_8h,f_8h,size); // f_8h <- (f_8h - A_8h * u_8h)
		T_8h.matrix_by_vector(f_8h, f_16h, size); // f_16h <- I_t_8h * f_8h
		iter = A_16h.gauss_seidel_criterion(e_16h,f_16h,criterion);  // Solve, A_16h * e_16h = f_16h , so we get e_16h
	
		// Interpolate the coarse-grid error to the fine grid //
	
		// 16h->8h //
		size = A_8h.r;
		I_8h.matrix_by_vector(e_16h, e_8h, A_16h.r);// e_8h <- I_8h * e_16h
		
		// 8h->4h //
		for(int i=0; i<size; i++)
			u_8h[i] += e_8h[i];						// u_8h <- u_8h + e_8h
		size = A_4h.r;
		I_4h.matrix_by_vector(u_8h, e_4h, A_8h.r);	// e_4h <- I_4h * e_8h
	
		// 4h->2h //
		for(int i=0; i<size; i++)
			u_4h[i] += e_4h[i];						// u_4h <- u_4h + e_4h
		size = A_2h.r;
		I_2h.matrix_by_vector(u_4h, e_2h, A_4h.r);	 // e_2h <- I_2h * e_4h
	
		// 2h->h //
		for(int i=0; i<size; i++)
			u_2h[i] += e_2h[i];						// u_2h <- u_2h + e_2h
		size = A_h.r;
		I_h.matrix_by_vector(u_2h, e_h, A_2h.r);	 // e_h <- I_h * u_2h
	
		// get u_h //
		for(int i=0; i<size; i++)
			u_h[i] += e_h[i];						// u_h <- u_h + e_h
	
		// Get error vector //
		for(int i=0; i<size; i++)
			error[i]=solution[i]-u_h[i];

		cout << " : Norm of error vector = "<<norm(error,size)<<endl;
	}
	
	// Calculate and Print the running time //
	end = clock();
	printf("\nTotal execution time : %.1lf (s)\n", (end-start)/(double)CLOCKS_PER_SEC);
	

	// Store the result vector u_h into the text file //
	size = A_h.r;
	ofstream outFile("../results/output.txt",n);
	for(int i=0;i<size-1;i++)
	{
		outFile << u_h[i] <<" ";
	}
	outFile << u_h[size-1];
	outFile.close();

	cout << endl;
	cout << "Draw the result plot using matlab function :" << endl;
	cout << "	drawResult()" << endl;
	cout << endl;
	cout << "Type anything to terminate" << endl;
	int Type_anything_to_terminate;
	cin >> Type_anything_to_terminate;
	return 0;
}
