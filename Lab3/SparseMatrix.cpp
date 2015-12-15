#include "stdafx.h"

SparseMatrix::SparseMatrix(int _rownum, int _colnum, int _nz) : rownum(_rownum), colnum(_colnum) {
	nz = _nz;
	vals = new double[nz];
	rows = new int[nz];
	cols = new int[nz];
}

SparseMatrix::~SparseMatrix() {
	delete[] vals;
	delete[] rows;
	delete[] cols;
}

void SparseMatrix::Generate(double value_lb, double value_hb) {
	for (int i = 0; i < nz; ++i) {
		vals[i] = (double)rand() / RAND_MAX * (value_hb - value_lb) + value_lb;
		rows[i] = rand() % rownum;
	}
	my_qsort(rows, 0, nz-1);

	int currowstart = 0;
	cols[0] = rand() % colnum;
	for (int ind = 1; ind < nz; ++ind) {
		if (rows[ind] != rows[currowstart]) {
			/* New row begins */
			my_qsort(cols, currowstart, ind-1);
			currowstart = ind;
			cols[ind] = rand() % colnum;
		}
		else {
			/* Generating unique column index */
			bool suitable = false;
			do {
				cols[ind] = rand() % colnum;
				suitable = true;
				for (int i = currowstart; i < ind; ++i)
					if (cols[i] == cols[ind]) {
						suitable = false;
						break;
					};
			}
			while (!suitable);
		}
	}
	my_qsort(cols, currowstart, nz-1);
}

void SparseMatrix::PrintAsArrays() {
	std::cout << "rows: ";
	for (int i = 0; i < nz; ++i)
		std::cout << rows[i] << " ";
	std::cout << "\n" << "cols: ";
	for (int i = 0; i < nz; ++i)
		std::cout << cols[i] << " ";
	std::cout << "\n" << "vals: ";
	for (int i = 0; i < nz; ++i)
		std::cout << vals[i] << " ";
	std::cout << "\n";
}

void SparseMatrix::SortRow(int row_start, int row_end) {
	my_qsort_bounded(cols, vals, row_start, row_end);
}

void SparseMatrix::ReallocateMemory(int actNZ) {
	int* _rows = new int[actNZ];
	int* _cols = new int[actNZ];
	double* _vals = new double[actNZ];
	for (int i = 0; i < (nz < actNZ ? nz : actNZ); ++i) {
		_rows[i] = rows[i];
		_cols[i] = cols[i];
		_vals[i] = vals[i];
	}
	delete[] rows;
	delete[] cols;
	delete[] vals;
	nz = actNZ;
	rows = _rows;
	cols = _cols;
	vals = _vals;
	//std::cout << "Memory reallocated!\n";
}

SparseMatrix& Multiply(SparseMatrix& A, SparseMatrix& B) {
	if (A.colnum != B.rownum) {
		throw;
	}
	int estNZ = A.nz*(B.nz/B.rownum+1);
	SparseMatrix& C = *(new SparseMatrix(A.rownum, B.colnum, estNZ));
	//std::cout << "estNZ = " << estNZ << "\n";
	int curNZ = 0;
	int currowstart = 0;
	for (int indA = 0; indA < A.nz; ++indA) {
		if (curNZ && C.rows[currowstart] != A.rows[indA]) {
			C.SortRow(currowstart, curNZ-1);
			currowstart = curNZ;
		}
		int indB = 0;
		while (indB < B.nz && B.rows[indB] < A.cols[indA])
			++indB;
		while (indB < B.nz && B.rows[indB] == A.cols[indA]) {
			int column = B.cols[indB];
			bool exists = false;
			for (int i = currowstart; i < curNZ; ++i)
				if (C.cols[i] == column) {
					exists = true;
					C.vals[i] += A.vals[indA]*B.vals[indB];
					//C.PrintAsArrays();
					break;
				}
			if (!exists) {
				if (curNZ == C.nz) C.ReallocateMemory(2*curNZ);
				C.rows[curNZ] = A.rows[indA];
				C.cols[curNZ] = column;
				C.vals[curNZ] = A.vals[indA]*B.vals[indB];
				++curNZ;
				//C.PrintAsArrays();
			}
			++indB;
		}
	}
	//std::cout << "actNZ = " << curNZ << "\n";
	C.SortRow(currowstart, curNZ-1);
	C.ReallocateMemory(curNZ);
	return C;
}

SparseMatrix& SparseMatrix::Scatter() {
	int ProcID, NumProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcID);
	MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);

	int q = nz / NumProcs;
	int r = nz % NumProcs;

	int* sendcounts = new int[NumProcs];
	int* displacements = new int[NumProcs];
	for (int k = 0; k < NumProcs - 1; ++k) {
		sendcounts[k] = q;
		displacements[k] = k*q;
	}
	sendcounts[NumProcs-1] = q+r;
	displacements[NumProcs-1] = (NumProcs-1)*q;

	SparseMatrix& part = *(new SparseMatrix(rownum, colnum, sendcounts[ProcID]));
	MPI_Scatterv(rows, sendcounts, displacements, MPI_INT, part.rows, sendcounts[ProcID], MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(cols, sendcounts, displacements, MPI_INT, part.cols, sendcounts[ProcID], MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(vals, sendcounts, displacements, MPI_DOUBLE, part.vals, sendcounts[ProcID], MPI_DOUBLE, 0, MPI_COMM_WORLD);

	delete[] sendcounts; delete[] displacements;
	return part;
}

void SparseMatrix::Broadcast() {
	int ProcID, NumProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcID);
	MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);

	MPI_Bcast(rows, nz, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(cols, nz, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(vals, nz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

SparseMatrix& SparseMatrix::Gather() {
	int ProcID, NumProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcID);
	MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);

	int estNZ = 0;
	int *recvcounts = 0, *displacements = 0;
	if (!ProcID) {
		recvcounts = new int[NumProcs];
		displacements = new int[NumProcs];
	}
	MPI_Gather(&nz, 1, MPI_INT, recvcounts+ProcID, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (!ProcID) 
		for (int k = 0; k < NumProcs; ++k) {
			displacements[k] = estNZ;
			estNZ += recvcounts[k];
		}

	int *_rows = 0, *_cols = 0;
	double *_vals = 0;
	if (!ProcID) {
		_rows = new int[estNZ];
		_cols = new int[estNZ];
		_vals = new double[estNZ];
	}
	MPI_Gatherv(rows, nz, MPI_INT, _rows, recvcounts, displacements, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gatherv(cols, nz, MPI_INT, _cols, recvcounts, displacements, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gatherv(vals, nz, MPI_DOUBLE, _vals, recvcounts, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	SparseMatrix& S = *(new SparseMatrix(rownum, colnum, estNZ));

	int curNZ = 0;
	int currowstart = 0;
	for (int ind = 0; ind < estNZ; ++ind) {
		if (_vals[ind]) {
			if (curNZ && _rows[ind] != S.rows[currowstart]) {
				S.SortRow(currowstart, curNZ-1);
				currowstart = curNZ;
			}
			bool exists = false;
			for (int i = currowstart; i < curNZ; ++i)
				if (S.cols[i] == _cols[ind]) {
						exists = true;
						S.vals[i] += _vals[ind];
						break;
				}
			if (!exists) {
					S.rows[curNZ] = _rows[ind];
					S.cols[curNZ] = _cols[ind];
					S.vals[curNZ] = _vals[ind];
					++curNZ;
			}
		}
	}
	S.ReallocateMemory(curNZ);
	if (!ProcID) {
		delete[] recvcounts; delete [] displacements;
		delete[] _rows; delete[] _cols; delete[] _vals;
	}
	return S;
}

bool AreEqual(SparseMatrix& A, SparseMatrix& B) {
	if (A.rownum == B.rownum &&
		A.colnum == B.colnum &&
		A.nz == B.nz) {
			bool res = true;
			for (int i = 0; i < A.nz; ++i)
				if (A.rows[i] != B.rows[i] ||
					A.cols[i] != B.cols[i] ||
					A.vals[i] != B.vals[i]) {
						res = false;
						break;
				}
			return res;
	}
	else return false;
}