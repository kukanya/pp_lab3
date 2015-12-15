#pragma once

class SparseMatrix {
private:
	const int rownum, colnum;
	int nz;
	int *vals, *rows, *cols;

	void SortRow(int row_start, int row_end);
	void ReallocateMemory(int actNZ);

public:
	SparseMatrix(int _rownum, int _colnum, int _nz);
	~SparseMatrix();
	void Generate();
	void PrintAsArrays();

	SparseMatrix& Scatter();
	void Broadcast();
	SparseMatrix& Gather();

	friend SparseMatrix& Multiply(SparseMatrix& A, SparseMatrix& B);
	friend bool AreEqual(SparseMatrix& A, SparseMatrix& B);
};