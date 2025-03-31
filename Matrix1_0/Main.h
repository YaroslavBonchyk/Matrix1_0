#pragma once

class Matrix {
private:
	int row, col;
	double** matrix;
public:

	//General
	Matrix(int row, int col);
	Matrix(const Matrix& M);
	~Matrix();
	void Fill();
	void FillRand();
	void Print();

	//Overload operators
	Matrix& operator = (const Matrix& A);
	friend Matrix operator +(const Matrix& A, const Matrix& B);
	friend Matrix operator -(const Matrix& A, const Matrix& B);
	friend Matrix operator *(const Matrix& A, const Matrix& B);
	friend Matrix operator *(double scalar, const Matrix& A);

	// Functions
	Matrix Transpose() const;
	double Determinant(int col) const;
	double Trace() const;
	double Minor(int col, int indrow, int indcol) const;
	Matrix Inverse(const Matrix& A);
	double DotProduct(const Matrix& V1, const Matrix& V2);
	double Norm()const;
	Matrix NormalizeVector(const Matrix& V);
	Matrix SplitVector(int w)const;
	Matrix SplitRow(int r)const;
	Matrix ZeroVector(int z)const;
	Matrix IdentityMatrix(int s) const;
	bool IsSpecial(const Matrix& A);

	//SLE
	Matrix REF(Matrix& RA);
	int Rank(Matrix& A);
	Matrix HSLE(Matrix& A);
	Matrix NHSLE(Matrix& A, Matrix& B);

	//QRD
	Matrix Qmatrix(const Matrix& A);
	Matrix Rmatrix(const Matrix& Q, const Matrix& A);
	Matrix QRalgBasic(const Matrix& A);
	Matrix HessenbergMatrix(const Matrix& A);
	Matrix QRalgHess(const Matrix& A);

	//EVD
	Matrix EigenValMatrix(const Matrix& A);
	Matrix EigenVectMatrix(const Matrix& A);

	//SVD
	Matrix Smatrix(const Matrix& A);
	Matrix Vmatrix(const Matrix& A);
	Matrix Umatrix(const Matrix& A, const Matrix& V, const Matrix& S);
};
