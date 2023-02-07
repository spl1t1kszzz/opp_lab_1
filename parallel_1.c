#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

const int n = 7;
const double t = 0.001;
const double epsilon = 1e-05;

void createMatrix(double *matrix, int size, double mainData, double otherData) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i == j)
                *(matrix + i * n + j) = mainData;
            else
                *(matrix + i * n + j) = otherData;
        }
    }
}

void printMatrix(const double *matrix, int size1, int size2) {
    for (int i = 0; i < size1; ++i) {
        for (int j = 0; j < size2; ++j) {
            printf("%lf ", *(matrix + i * n + j));
        }
        printf("\n");
    }
}

double norm(const double *vector, int size) {
    double norm = 0.0;
    for (int i = 0; i < size; ++i) {
        norm += *(vector + i) * *(vector + i);
    }
    return sqrt(norm);
}

int main(int argc, char **argv) {
    int rank, size;
    double *matrix = NULL;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (rank == 0) {
        matrix = calloc(n * n, sizeof(double));
        createMatrix(matrix, n, 2.0, 1.0);
        for (int i = 1; i < size; ++i) {
            MPI_Send(matrix + n * (n / size + n % size) + (n/size) * n * (i - 1), n * (n / size), MPI_DOUBLE, i, 123, MPI_COMM_WORLD);
        }
        printMatrix(matrix, (n / size + n % size), n);
    }
    else {
        matrix = calloc(n * (n / size), sizeof(double));;
        if (MPI_Recv(matrix , n * (n / size), MPI_DOUBLE, MPI_ANY_SOURCE, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE) == MPI_SUCCESS) {
            printMatrix(matrix, n / size , n);
        }
        else {
            printf("Error!");
        }
    }
    MPI_Finalize();
    return 0;
}
