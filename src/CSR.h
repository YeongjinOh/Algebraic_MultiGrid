#include <iostream>

using namespace std;

class CSR {
//private:
public:
	int size; // Size of array. it should be greater than the number of non-zero elements.
	int cnt;  // Count how many elements have been inserted. 

	int r; // size of row of matrix
	int c; // size of column of matrix

	int *row; // The array that contains the index of row
	int *col; // The array that contains the index of column
	double *data; // The array that contains the data of matrix.
	// For example row[0]=5, col[0]=3, data[0]=0.1 means, the (5,3)-th value of matrix is 0.1
	
//public:
	CSR();
	~CSR();
	CSR(const CSR& oldCSR);

	/* when we build CSR matrix, we need only these two function */
	CSR(int size_);
	CSR(int size_, int r_, int c_);
	void SetValue(int r, int c, double value);

	void    Show();
	void	Show_ith_value(int i);
	void	Info();
	void	print();
	void	print(int part_r, int part_c);

	void    gauss_seidel(double *u, double *f, int n) ;
	int     gauss_seidel_criterion(double * u, double * f, double criterion);
	void	build_I(CSR A, int *c, int *d, int *index, int n);
	void    matrix_by_vector(double *u, double *f, int size);
	void    transpose();
	void    copy_CSR(CSR *copy);
	void    multplication(CSR *src1, CSR *src2);

	double max();
	double min();
	double row_minimum(int i);
	double ij_value(int i, int j);
};

CSR::CSR() {
    size   = 0;
	cnt    = 0;
	row    = 0;
    col    = 0;
    data   = 0;
	r = 0;
	c = 0;
}
CSR::~CSR() {
    row    = 0;
	col    = 0;
	data    = 0;
}
CSR::CSR(const CSR& oldCSR)
    {
		size = oldCSR.size;
		cnt = oldCSR.cnt;
		row = new int[size];
        col = new int[size];
		data = new double[size];
		for(int i=0; i<size; i++)
		{
				row[i] = oldCSR.row[i];
				col[i] = oldCSR.col[i];
				data[i] = oldCSR.data[i];
		}
    }
CSR::CSR(int size_) {
    size = size_;
	cnt = 0; // Initialize cnt as 0.
	row = new int[size_];
    col = new int[size_];
    data = new double[size_];
	/* Since we know how many elements are available (by cnt), we don't access to and care about the remaining garbage values. */
	for(int i=0; i<size; i++)
	{
		row[i] = 0;
		col[i] = 0;
		data[i] = 0;
	} 
}
CSR::CSR(int size_, int r_, int c_){
    size = size_;
	r=r_;
	c=c_;
	cnt = 0; // Initialize cnt as 0.
	row = new int[size_];
    col = new int[size_];
    data = new double[size_];
	/* Since we know how many elements are available (by cnt), we don't access to and care about the remaining garbage values. */
	for(int i=0; i<size; i++)
	{
		row[i] = 0;
		col[i] = 0;
		data[i] = 0;
	} 
}
void CSR::SetValue(int r, int c, double value) // Set the (r,c)-th value of matrix as 'value'
{
	if(cnt>=size) // exception
	{
		std::cout << "Type anything to terminate" << std::endl; 
	}
	else
	{
		row[cnt] = r;
		col[cnt] = c;
		data[cnt] = value;
		
		cnt++;
	}
}
void CSR::Show() { // Show all the elements in array row, col,and data
	
	printf("row :");
	for(int i=0; i<cnt; i++) 
		printf(" %d", row[i]);
	printf("\n");

	printf("col :");
	for(int i=0; i<cnt; i++) 
		printf(" %d", col[i]);
	printf("\n");

	printf("values :");
	for(int i=0; i<cnt; i++) 
		printf(" %.2lf", data[i]);
	printf("\n");

}
void CSR::Show_ith_value(int i) { // Show the information of i-th inserted value
	
	printf("row = %d \n",row[i-1]);
	printf("col = %d \n",col[i-1]);
	printf("val = %.1lf \n",data[i-1]);
}
double CSR::max() // Return the maximum value of the matrix
{
	if(cnt == 0)
	{
		cout << "** The matrix doesn't have any element! **" <<endl;
		return 0; // exception
	}
	else
	{
		double max = data[0];
		for(int i=1; i<cnt; i++)
			if(max<data[i])
				max=data[i];

		return max;

	}
}
double CSR::min() // Return the minimum value of the matrix
{
	if(cnt == 0)
	{
		cout << "** The matrix doesn't have any element! **" <<endl;
		return 0; // exception
	}
	else
	{
		double min = data[0];
		for(int i=1; i<cnt; i++)
			if(min>data[i])
				min=data[i];

		return min;

	}
}
double CSR::row_minimum(int i) // Find the minimum value of i-th row
{
	// Assume that there exist at least one element less than 100 //
	double min=100;
	int cur=0;
	while(row[cur] < i)
		cur++;
	if(row[cur] == i)
	{
		while(row[cur] == i)
		{
			if(min>data[cur])
				min=data[cur];
			cur++;
		}
	} else  // i-th row doesn't have any element.
	{
		cout << "** " << i << "-th row doesn't have any element **" << endl;
	}
	return min;
}
double CSR::ij_value(int i, int j) // Return the (i,j)-th value
{
	for(int cur=0; cur<cnt; cur++)
	{
		if(row[cur]==i && col[cur]==j)
		{
			return data[cur];
		} else {
				cout << "** The matrix doesn't have non-zero element at ("<< i << ", " << j << ") **" <<endl;
		}
	}
	return 0;
}
void CSR::Info() // Show the information about size
{
	cout << "r = "<<r<<", c = "<<c<<", cnt = "<<cnt<<", size = "<<size<<endl;	
	cout << endl;
}
void CSR::print() // Print the matrix
{
	int cur=0;
	for(int i=1; i<=r; i++)
	{
		for(int j=1; j<=c; j++)
		{
				if(row[cur] == i && col[cur] == j)
				{
					if(data[cur]<0)
						printf("%.1lf ", data[cur]);
					else
						printf(" %.1lf ", data[cur]);
					cur++;
				} else {
					cout<< "  0  ";
				}
		}
	cout << endl;
	}
}
void CSR::print(int part_r, int part_c)	// Print the upper-left part_r by part_c submatrix of this matrix
{
	int cur=0;
	for(int i=1; i<=part_r; i++)
	{
		while(row[cur] < i)
			cur++;
		for(int j=1; j<=part_c; j++)
		{
				if(row[cur] == i && col[cur] == j)
				{
					if(data[cur]<0)
						printf("%.2lf ", data[cur]);
					else
						printf(" %.2lf ", data[cur]);
					cur++;
				} else {
					cout<< "  0   ";
				}
		}
	cout << endl;
	}
}
void CSR::gauss_seidel(double * u, double * f, int n) // Solve CSR*u=f by using gauss_seidel n times.
{
	for(int k=0; k<n; k++)
	{
		int cur=0;
		double new_sum, old_sum;
		double a_ii;
		for(int i=1; i<=r; i++) // r = length of row of A = n*n
		{
			new_sum=0;
			old_sum=0;
			a_ii=0;

			while(row[cur]<i) // Make cur point to row i
				cur++;

			while(row[cur] == i)
			{
				if(col[cur] < i)
				{
					new_sum += data[cur]*u[col[cur]-1];
				} else if(col[cur] > i)
				{
					old_sum += data[cur]*u[col[cur]-1];
				} else {
					a_ii = data[cur];
				}
				cur++;
			}
			if(a_ii == 0) // Check whether A is positive diagonal matrix.
			{
				for(int m=0; m<100; m++)
					cout << "The matrix has zero element on the diagonal!";
				cout << "Check in Gauss_seidel in Multigrid.h" <<endl;
			}
			u[i-1] = (f[i-1] - new_sum - old_sum)/a_ii;
		}
	}
}
int CSR::gauss_seidel_criterion(double * u, double * f, double criterion)
{
	int iter=0;
	double norm_value = 100;
	double norm_divisor = 0;
	double *old_u = new double[r];
	double *norm_u = new double[r];
	double new_sum,old_sum, a_ii;
			

	for(int i=0; i<r; i++)
		old_u[i]=u[i];

	while(norm_value > criterion)
	{
		int cur=0;
		for(int i=1; i<=r; i++) // r = length of row of A = n*n = size of u 
		{
			new_sum=0;
			old_sum=0;
			a_ii=0;

			while(row[cur]<i)
				cur++;
			while(row[cur] == i)
			{
				if(col[cur] < i)
				{
					new_sum += data[cur]*u[col[cur]-1];
				} else if(col[cur] > i)
				{
					old_sum += data[cur]*u[col[cur]-1];
				} else {
					a_ii = data[cur];
				}
				cur++;
			}
			u[i-1] = (f[i-1] - new_sum - old_sum)/a_ii;
		}
		for(int i=0; i<r; i++)
		{
			norm_u[i]=u[i]-old_u[i];
			old_u[i]=u[i];
		}
		
		norm_value=0;
		for(int i=0; i<r; i++)
		{
			norm_value += (norm_u[i]*norm_u[i]);
		}

		norm_divisor = 0;
		for(int i=0; i<r; i++)
		{
			norm_divisor += u[i]*u[i];
		}
		//norm_value /= (double)r; // norm_value = || old_u - u || 
		 norm_value /= norm_divisor; // norm_value = || old_u - u || / || u ||
		iter++;
	}
	//cout << "The last norm_value = "<<norm_value <<endl;
	return iter;
}
void CSR::build_I(CSR A, int * c, int * d, int * index, int size) // Make this matrix interpolation operator I by using matrix A, array c and d
{
	int c_cnt = 1;
	int * c_idx = new int[size]; // Where C will be located
	int idx, divisor;

	for(int i=0; i<size; i++)
	{
		if (c[i] == 1)
		{
			c_idx[i] = c_cnt;
			c_cnt++;
		}
	}

	c_cnt = 1;
	for(int i=0; i<size; i++)
	{
		if(c[i] == 1) // In the case of C
		{
			row[cnt] = i+1;
			col[cnt] = c_cnt;
			data[cnt] = 1;
		
			cnt++;
			c_cnt++;
		}
		else // Case of F
		{
			divisor = 0; // divisor counts the number of C dependent to F
			for(int k=0; k<d[i]; k++)
			{
				idx = index[k*size+i]; // index of neighbor of F
				if(c[idx] == 1)
					divisor++;

			}
			for(int k=0; k<d[i]; k++)
			{
				idx = index[k*size+i]; // index of neighbor of F
				if(c[idx] == 1)
				{
					row[cnt] = i+1;
					col[cnt] =  c_idx[idx];
					data[cnt] = 1/(double)divisor; // It might need to be modified.
					cnt++;
				}
			}
		}
	}
}
void CSR::matrix_by_vector(double * u, double * f, int size) // Assign CSR*u into f, size = length of u
{
	if(size != c)
	{
		cout << "Wrong input size in matrix_by_vector!!" << endl;
		cout << "r = " << r <<", size = " << size << endl;
		return;
	}
	int cur = 0;
	for(int i=1; i<=r; i++)
	{
		while(row[cur]<i)
				cur++;
		f[i-1]=0;
		while(row[cur] == i)
		{
			//cout << "f[i], i="<<i-1<< "// u[j], j="<<col[cur]-1<<endl; 
			f[i-1] += data[cur]*u[col[cur]-1];
			cur++;
		}


	}

}
void CSR::copy_CSR(CSR * copy) // Get copied from 'copy'
{
	size=copy->size;
	cnt=copy->cnt;
	r=copy->r;
	c=copy->c;
	for(int i=0; i<cnt; i++)
	{
		row[i]=copy->row[i];
		col[i]=copy->col[i];
		data[i]=copy->data[i];
	}
}
void CSR::transpose() // Transpose this matrix
{
	int * tmp_array;
	// Swap row and column //
	tmp_array=row;
	row=col;
	col=tmp_array;

	int tmp;
	tmp = r;
	r = c;
	c = tmp;
	double tmp_d;
	// Bubble sorting //
	for(int i=cnt-2; i>=0; i--) // Assume that cnt>=2
	{
		for(int j=i; j<cnt-1; j++)
		{
			if(row[j] > row[j+1])
			{
				tmp=row[j+1];
				row[j+1]=row[j];
				row[j]=tmp;
				tmp=col[j+1];
				col[j+1]=col[j];
				col[j]=tmp;
				tmp_d=data[j+1];
				data[j+1]=data[j];
				data[j]=tmp_d;
			}
		}
	}
}
void CSR::multplication(CSR *src1, CSR *src2) // Assign src1*transpose(src2) into this matrix.
	// Note that src2 is transposed form of right matrix of multiplication in order to use CSR property easily.
{
	if(src1->c != src2->c)
	{
		cout<< "Wrong input size in matrix multiplication!!" <<endl;
	} else
	{
		r = src1->r;
		c = src2->r;
		cnt=0;
		double sum;
		int cur_i=0, cur_j=0, store_i;
		
		for (int i=1; i<=r; i++)
		{
			while(src1->row[cur_i]<i)
				cur_i++;
			store_i=cur_i;
			for(int j=1; j<=c; j++)
			{
				cur_i = store_i;
				sum=0;
				cur_j=0;
				while(src2->row[cur_j]<j)
				{
				//	cout << "row = "<<row[cur_j] << ", j = "<<j<<", cur_j = "<<cur_j<<endl;				
					cur_j++;
				}
				while(src1->row[cur_i] == i && src2->row[cur_j] == j)
				{
				//	cout << "i,j = " << i<< ", "<<j<< "  cur : "<< cur_i << " " << cur_j<<endl;
					if(src1->col[cur_i]>src2->col[cur_j])
					{
						cur_j++;
					} else if(src1->col[cur_i]<src2->col[cur_j]) {
						cur_i++;
					} else	{
					//	cout << "i,j = " << i<< ", "<<j<< "  cur : "<< cur_i << " " << cur_j<<endl;
						sum+= src1->data[cur_i]*src2->data[cur_j];
						cur_i++;
						cur_j++;					
					}
					//cout << "i,j = " << i<< ", "<<j<< "  cur : "<< cur_i << " " << cur_j<<endl;
				}
				if(sum != 0)
					SetValue(i,j,sum);
				
			}
		}
	}
}