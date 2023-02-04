#include "mpi.h"
#include <stdio.h>

int main(int argc, char** argv) {
    int size, rank;
    double x[3] = {1.2, 5.5, 6.8};
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (size >= 3) {
        if (rank == 0) {
            MPI_Send(x, 3, MPI_DOUBLE, 2, 123, MPI_COMM_WORLD);
        }
        if (rank == 2) {
            if (MPI_Recv(x, 3, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE) == MPI_SUCCESS) {
                for (int i = 0; i < 3; ++i) {
                    printf("%lf ", x[i]);
                }
        }}
    }
    MPI_Finalize();
    return 0;


}