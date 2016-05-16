#include "EigenWrapper.h"
#include <iostream>
using namespace std;

#define INITIAL_N 128;	// The length of row and column of Initial block.
#define THRESHOLD 0.01;	// The threashold value of strongly dependency. It must be greater than 0 and less than or equal to 1
#define USE_THRESHOLD false;
#define ITERATION 3;	// The number of iteration in Gauss_seidel
#define MAX_D 20;		// The maximal number of dependancies
#define CYCLE 10;		// V-cycle iteration number.
#define TRACKING true;	// If you want to check what's being calculated, Set it true.



void print_i(int n, int *d) // Print int array of nxn matrix 
{
	for(int i=0; i<n*n; i++)
	{
		cout << d[i] << " ";
		if(i%n == n-1)
			cout << endl;
	}
}
void print_d(int n, double *d) // Print double array of nxn matrix //
{
	for(int i=0; i<n*n; i++)
	{
		cout << d[i] << " ";
		if(i%n == n-1)
			cout << endl;
	}
}
void print_d_mk(int n, double *d, int m, int k) // print the first m by k submatrix of d
{
	for(int i=1; i<=m; i++)
	{
		for(int j=1; j<=k; j++)
			cout << d[n*(i-1) + (j-1)] << " ";
		cout << endl;
	}
}
void print_vector(int size, int *v) // Print the first 'size' elements of int vector v 
{
	for(int i=0; i<size; i++)
		cout<<" "<<v[i]<<endl;
}
void print_vector_d(int size, double *v) // Print the first 'size' elements of double vector v 
{
	for(int i=0; i<size; i++)
		cout<<" "<<v[i]<<endl;
}
double max(double * arr, int size)	// Get maximum value of given double array(or vector) with length of 'size'
{
	double max = arr[0];
	for(int i=1; i<size; i++)
	{
		if(max<arr[i])
			max=arr[i];
	}

	return max;
}
double average(double * arr, int size) // Get average value of given double array(or vector) with length of 'size'
{
	double avg = 0;
	for(int i=0; i<size; i++)
	{
		avg += arr[i];
	}

	return avg/size;
}
double norm(double * arr, int size) // Get euclidean norm of given double array(or vector) with length of 'size'
{
	double norm = 0;
	for(int i=0; i<size; i++)
	{
		norm += arr[i]*arr[i];
	}

	return norm/size;
}

void copy_int_array(int * a, int * b, int size) // Copy int array a into array b
{
	for(int i=0; i<size; i++)
		b[i]=a[i];
}
void make_f_remainder(CSR * A, double * u, double * f, int size) // f <- (f-Au), size = length of u
{
	bool print = false;

		double * sub = new double[size];
		A->matrix_by_vector(u,sub,size); // sub <- A*u
	if(print)
	{
		A->print();
		cout<<"	  *   "<<endl;
		print_vector_d(size,u);
		cout<<"   =   "<<endl;
		print_vector_d(size,sub);
		cout << " f = " <<endl;
		print_vector_d(size,f);
	}
		for(int i=0; i<size; i++)
			f[i] -= sub[i];
	if(print)
	{
		cout << " r = " <<endl;
		print_vector_d(size,f);
	}
}

int find_first_zero(int * c, int k) // Find the first index of array c whose value is zero
{
	for(int i=0; i<k; i++)
		if(c[i]==0)
			return i;
	return -1;
}
int maximum_index(int *d, int *old_d, int *c, int *index, int size, int C_id)
{
	// At the first time, we search all of the dependency values and find the maximum //
	// After that, we search the maximum dependency among the elements that are strongly dependent to F //
	// If all of them are already assigned with C value, we search all of the dependency values and find the maximum //

	int max;
	int idx=-1;
	
	if( C_id == -1)
	{
		idx = find_first_zero(c,size);
		max = d[idx];
		for(int i=idx; i<size; i++)	
		{
			if(max < d[i] && c[i]==0)
			{
				max = d[i];
				idx = i;
			}
		}
	} else {
		max = -1;
		for(int i=0; i<old_d[C_id]; i++)
		{
			int F_id = index[i*size+C_id]; // The index of the candidate who would be F. (or already F)
			for(int j=0; j<old_d[F_id]; j++)
			{
				int id = index[j*size+F_id]; // The index of the element dependent to F.
				if(max < d[id] && c[id]==0)
				{
					max = d[id];
					idx = id;
				}
			}
		}
		if(idx < 0)
		{
			idx = find_first_zero(c,size);
			max = d[idx];
			for(int i=idx; i<size; i++)	
			{
				if(max < d[i] && c[i]==0)
				{
					max = d[i];
					idx = i;
				}
			}
		}
	}

	if(idx < 0) // Exception handling for the case never happen
		throw runtime_error(" *** Negative index in maximum_index function!! *** ");

	return idx;
}
void increase_neighbor(int *d, int n, int idx) // idx is index of F, So we will increase the value d of neighbors of F
{
	//     F is    //
	if(idx%n != 0) // Not in the first column
		d[idx-1]++;
	if(idx%n != n-1) // Not in the last column
		d[idx+1]++;
	if(idx >= n) // Not in the first row
		d[idx-n]++;
	if(idx < n*n-n) // Not in the last row
		d[idx+n]++;
}
void build_c_manually(CSR csr, int * d, int * c, int n) // We build initial c manually. (size of c and d are n*n)
{
	for(int i=1; i<=n; i++)
		for(int j=1; j<=n; j++)
			c[(i-1)*n + (j-1)] = (i+j)%2+1;
}
int build_c(CSR csr, int *d, int *old_d, int *c, int *index, int size) // Bulid Coarsing array c
{
	// The value of c means		//
	// if 0, it's not-assigned	//
	// if 1, then C				//
	// if 2, then F				//

	int prev_idx = -1;
	
	// Initialize array c //
	for(int i=0; i<size; i++)
	{
		c[i]=0;
	}
try {
	int cnt =0;
	int C_cnt = 0; // Count the number of C's	
	while(cnt < size)
	{
		
			int C_id = maximum_index(d, old_d, c, index, size, prev_idx);
		
		prev_idx = C_id;
		c[C_id] = 1;
		cnt++; // Counts the number of C's and F's
		C_cnt++;
		for(int i=0; i<old_d[C_id]; i++)
		{
			int F_id = index[i*size+C_id]; // The index of the candidate who would be F. (or already F)
			if( c[F_id] == 0)
			{
				c[F_id] = 2; // Set it F
				cnt++;
				for(int j=0; j<old_d[F_id]; j++)
				{
					int id = index[j*size+F_id]; // The index of the candidate who would be F. (or already F)
					d[id]++;
				}
			}
		}
	}
		return C_cnt;
	} catch (exception const& e) {
		cout << "\nException: " << e.what() << "\n";

	}

}

void build_d (CSR csr,int size, int * d, int * index) // Build dependency array d by using matrix A
{
	int current_cnt; // Count how many elements are dependent in this row.
	int row_array_index=0;
	double threshold = THRESHOLD;
	bool use_threshold = USE_THRESHOLD;
	int max_D = MAX_D;
	for (int row=1; row<=size; row++)
	{
		current_cnt = 0;
		while(csr.row[row_array_index]==row)
		{
			if(use_threshold)
			{
				if(csr.col[row_array_index] != row && csr.data[row_array_index] <= threshold*csr.row_minimum(row)) 
				// It takes so long time, so we set the threshold value small enough so that we can get the same result without getting minimum value of row //
				{
					index[current_cnt*size+ (row-1)] = csr.col[row_array_index]-1 ;	// Since column starts from 1. we should minus 1
					current_cnt++;
				}
			} else {
				if(csr.col[row_array_index] != row && csr.data[row_array_index] < 0) 
				{
					index[current_cnt*size+ (row-1)] = csr.col[row_array_index]-1 ;	// Since column starts from 1. we should minus 1
					current_cnt++;
				}
			}
			row_array_index++;
		}
		if(current_cnt>max_D)
			cout<<" ** dependency value("<<current_cnt<<") is greater than max_D("<<max_D<<")! **" <<endl;
		d[row-1]=current_cnt;
	}
}
void build_A(CSR * csr, int n) // Build the initial n*n by n*n matrix A
{
	
	// In the first block row  [ D I O O O ... O ]
	csr->SetValue(1,1,4);
	csr->SetValue(1,2,-1);
	csr->SetValue(1,n+1,-1);
	for(int j=2; j<=n-1; j++) // In j-th row of the first block row, we insert 4 elements in each row. 
	{
		csr->SetValue(j,j-1,-1);
		csr->SetValue(j,j,4);
		csr->SetValue(j,j+1,-1);
		csr->SetValue(j,j+n,-1);
	}
	csr->SetValue(n,n-1,-1);
	csr->SetValue(n,n,4);
	csr->SetValue(n,n+n,-1);

	/* In the middle block rows [ I D I O ... O]
								[ O I D I ... O]
								 ...
								[ O ...   I D I]
	*/
	for(int i=2; i<=n-1; i++) // int i-th block row
	{
		csr->SetValue((i-1)*n+1,(i-2)*n+1,-1);
		csr->SetValue((i-1)*n+1,(i-1)*n+1,4);
		csr->SetValue((i-1)*n+1,(i-1)*n+2,-1);
		csr->SetValue((i-1)*n+1,	i*n+1,-1);
		for(int j=2; j<=n-1; j++) // In j-th row of the first block row, we insert 4 elements in each row. 
		{
			csr->SetValue((i-1)*n+j,(i-2)*n+j,	-1);
			csr->SetValue((i-1)*n+j,(i-1)*n+j-1,-1);
			csr->SetValue((i-1)*n+j,(i-1)*n+j	,4);
			csr->SetValue((i-1)*n+j,(i-1)*n+j+1,-1);
			csr->SetValue((i-1)*n+j,	i*n+j	,-1);
		}
		csr->SetValue((i-1)*n+n,(i-2)*n+n,	-1);
		csr->SetValue((i-1)*n+n,(i-1)*n+n-1,-1);
		csr->SetValue((i-1)*n+n,(i-1)*n+n,	4);
		csr->SetValue((i-1)*n+n,	i*n+n,	-1);
	}
		
	
	// In the first block row  [ D I O O O ... O ]
	csr->SetValue((n-1)*n+1,(n-2)*n+1,-1);
	csr->SetValue((n-1)*n+1,(n-1)*n+1,4);
	csr->SetValue((n-1)*n+1,(n-1)*n+2,-1);
	for(int j=2; j<=n-1; j++) // In j-th row of the first block row, we insert 4 elements in each row. 
	{
		csr->SetValue((n-1)*n+j,(n-2)*n+j,	-1);
		csr->SetValue((n-1)*n+j,(n-1)*n+j-1,-1);
		csr->SetValue((n-1)*n+j,(n-1)*n+j	,4);
		csr->SetValue((n-1)*n+j,(n-1)*n+j+1,-1);
	}
	csr->SetValue((n-1)*n+n,(n-2)*n+n,	-1);
	csr->SetValue((n-1)*n+n,(n-1)*n+n-1,-1);
	csr->SetValue((n-1)*n+n,(n-1)*n+n,	4);

}
void build_I(CSR * I, CSR A, int * c, int n) // Build Interpolation operator I using matrix A and array c
{
	int cnt = 1;
	int * index = new int[n*n]; // where C will be located
	for(int i=0; i<n*n; i++)
	{
		if (c[i]==1)
		{
			index[i]=cnt;
			cnt++;
		}
	}
	cnt = 1;
	for(int i=0; i<n*n; i++)
	{
		if(c[i] == 1) // In the case of C
		{
			I->SetValue(i+1,cnt,1);
			cnt++;
		}
		else // Case of F
		{
			for(int k=0; k<n*n; k++)
			{
				if(A.ij_value(i+1 , k+1) < 0)
				{
					I->SetValue(i+1 , index[k] ,0.25);
				}
				
 
			}
		}
	}

}
void build_f(double * f, int n) // Calculate the vector b(h*h*f) and store it into the array f
{
	double delta = 1/(double)(n+1);
	for(int r=1; r<=n; r++) 
		for(int c=1; c<=n; c++)
		{
			double y = (double)r*delta;
			double x = (double)c*delta;
			f[(r-1)*n + (c-1)] = 2*((1-6*x*x)* y*y *(1-y*y) + (1-6*y*y)* x*x *(1-x*x))*delta*delta  ; //(r,c)-th element
		}
}
void build_solution(double * solution, int n) // Calculate the exact solution and store it into the array
{
	for(int r=1; r<=n; r++) 
		for(int c=1; c<=n; c++)
		{
			double y = (double)r/(double)(n+1);
			double x = (double)c/(double)(n+1);
			solution[(r-1)*n + (c-1)] = (x*x - x*x*x*x)*(y*y*y*y - y*y); //(r,c)-th element
		}
}
int coarse(/* Inputs: */CSR * A, int size, /* outputs: */ CSR * I, CSR * T) // Build interpolation operator I and T by using A
{
	int * d = new int[size];
	int * old_d = new int[size];
	int * c = new int[size];
	int max_D = MAX_D; // Maximal number of non-zero elements in row
	int * index = new int[max_D*size]; // the first n*n contains the first dependent element of u
	// For example, if u[5] is strongly dependent with u[1], u[3], u[7], u[9], index[5]=1, index[n*n+5]=3, index[2*n*n+5]=7, index[3*n*n+5]=9
	
	bool tracking = TRACKING;
	
	build_d(*A,size,d,index);
	copy_int_array(d,old_d,size); // Should store d before being changed //
	

	if(tracking)
		cout << "Building C...";
	int C_cnt;
	int n = INITIAL_N;
	
	if (size == n*n) // To save time, use manually calculated C only at the first time.
	{
		build_c_manually(*A,d,c,n);
		C_cnt = size/2;
	} else {
		C_cnt = build_c(*A,d,old_d,c,index,size);		 // Takes O(size^2)
	}
	if(tracking)
		cout << "Done / ";

	int non_zeros = max_D*size;
	*I = CSR(non_zeros,size,C_cnt); // I is 512*512 by 512*256 matrix
	*T = CSR(non_zeros,C_cnt,size);
	
	if(tracking)
		cout << "Building I...";
	I->build_I(*A, c, old_d, index, size);	 // Takes O(size)
	if(tracking)
		cout << "Done / ";
	
	transpose(I,T); // T <- trans(I)
	
	return C_cnt;
}
void restrict(/* Inputs: */CSR * A, double * f, int size, /* outputs: */ double * u, double * f_2h, CSR * T)
{
	int iteration = ITERATION;

	bool print = false;


	// Relax 'iteration' times on Au=f with initial guess 'u' //
	A->gauss_seidel(u,f,iteration); // u <- G(A,f)
	make_f_remainder(A,u,f,size); // f <- (f-Au)
	
	if(print)
	{	
		cout<< " Remainder r=f-Au "<< endl;
		print_d_mk(1,f,10,1);
		cout<<endl;
	}

	T->matrix_by_vector(f, f_2h, size); // f_2h <- I_t*f
}
void coarse_A(CSR *T, CSR *A, CSR *A_2h) // use CSR.h, so slow
{
	CSR temp=CSR(A->size);
	temp.multplication(T,A);
	A_2h->multplication(&temp,T);
} 
void coarse_A(CSR *T, CSR *A, CSR *I, CSR *A_2h)
{
	bool tracking = TRACKING;
	if(tracking)
		cout << "Building A...";

	CSR temp=CSR(A->size);
	csr_mult(T,A,&temp);
	csr_mult(&temp,I,A_2h);
	
	if(tracking)
		cout << "Done" << endl;
}
