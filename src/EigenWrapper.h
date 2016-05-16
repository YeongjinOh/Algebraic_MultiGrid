#include <Eigen\SparseCore>
#include <vector>
#include <iostream>
#include "CSR.h"


typedef Eigen::SparseMatrix<double, 1> SpMat;
typedef Eigen::Triplet<double> triple;

SpMat our2lib_csr_obsolute(CSR *src)
{
	std::vector<triple> coefficients;

	for (int i = 0; i < src->size; i++)
	{
		if (src->row[i] == 0 || src->col[i] == 0)
			continue;
		triple newCoff(src->row[i]-1, src->col[i]-1, src->data[i]);
		coefficients.push_back(newCoff);
	}
	
	/*
	for (int i = 0; i < coefficients.size(); i++)
	{
		std::cout << coefficients[i].row() << ", " << coefficients[i].col() << ", " << coefficients[i].value() << endl;
	}
	*/

	SpMat res(src->r, src->c);
	res.setFromTriplets(coefficients.begin(), coefficients.end());
	return res;
}

SpMat our2lib_csr(CSR *src)
{
	SpMat res(src->r, src->c);
	for (int i = 0; i < src->size; i++)
	{
		res.insert(src->row[i] - 1, src->col[i] - 1) = src->data[i];
	}
	return res;
}

CSR lib2our_csr(SpMat *src)
{
	CSR res(src->nonZeros(), src->rows(), src->cols());
	for (int k = 0; k < src->outerSize(); k++)
	{
		for (SpMat::InnerIterator it(*src, k); it; it++)
		{
			res.SetValue(it.row()+1, it.col()+1, it.value());
		}
	}
	return res;
}

void csr_mult(CSR *src1, CSR *src2, CSR *res)
{

	//cout << "our2lib..."<< endl;
	SpMat mat1 = our2lib_csr_obsolute(src1);
	SpMat mat2 = our2lib_csr_obsolute(src2);
	//cout << "Multiplication..."<< endl;
	SpMat mat_mult = mat1 * mat2;
	//cout << "lib2our..."<< endl;
	*res = lib2our_csr(&mat_mult);
	//cout << "Done"<< endl << endl;
	res->c=src2->c;
	res->r=src1->r;
}



void transpose(CSR *src, CSR *res)
{
	SpMat mat = our2lib_csr_obsolute(src);
	SpMat trans = mat.transpose();
	*res = lib2our_csr(&trans);
	res->c=src->r;
	res->r=src->c;
}

