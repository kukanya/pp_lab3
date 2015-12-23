// Lab3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

int main(int argc, char* argv[])
{
	const int m = atoi(argv[1]);
	const int nz = atoi(argv[2]);
	const int ExpNumber = (argc > 3) ? atoi(argv[3]) : 1;
	const double lb = (argc == 6) ? atof(argv[4]) : 1.0;
	const double hb = (argc == 6) ? atof(argv[5]) : 9.0;

	bool tiny = (m <= 7 ? true : false);

	double serialStartTime, serialTotalTime, 
		   parallelStartTime, parallelTotalTime;
	srand(time(0));

	int InitRes, ProcID;
	if (InitRes = MPI_Init(&argc, &argv)) {
		std::cout << "Error when initialising MPI, return code = " << InitRes << "\n";
		MPI_Abort(MPI_COMM_WORLD, InitRes);
	}
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcID);

	for (int i = 0; i<ExpNumber; ++i) {

		SparseMatrix A(m, m, nz);
		SparseMatrix B(m, m, nz);

		if (!ProcID) {
			std::cout << "=== Experiment " << i+1 << " ===\n";

			A.Generate(lb, hb);
			if (tiny) {
				std::cout << "matrix A:\n";
				A.PrintAsArrays();
				std::cout << "\n";
			}
			B.Generate(lb, hb);
			if (tiny) {
				std::cout << "matrix B:\n";
				B.PrintAsArrays();
				std::cout << "\n";
			}
		}
	
		parallelStartTime = MPI_Wtime();

		SparseMatrix partA = A.Scatter();
		B.Broadcast();
		SparseMatrix partC = Multiply(partA, B);
		SparseMatrix ParallelRes = partC.Gather();

		parallelTotalTime = MPI_Wtime() - parallelStartTime;

		if (!ProcID) {

			serialStartTime = MPI_Wtime();
			SparseMatrix SerialRes = Multiply(A, B);
			serialTotalTime = MPI_Wtime() - serialStartTime;

			if (AreEqual(SerialRes, ParallelRes))
				std::cout << "Results are equal\n";
			else 
				std:: cout << "Results are NOT equal!\n";
			if (tiny) {
				std::cout << "matrix C:\n";
				SerialRes.PrintAsArrays();
				std::cout << "\n";
			}
			else 
				std::cout << "Serial:   " << serialTotalTime << "s\n";
			if (tiny) {
				std::cout << "matrix PC:\n";
				ParallelRes.PrintAsArrays();
				std::cout << "\n";
			}
			else
				std::cout << "Parallel: " << parallelTotalTime << "s\n\n";
		}
	}
	MPI_Finalize();
	return 0;
}

