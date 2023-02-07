#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

const int n = 7;
const double t = 0.001;
const double epsilon = 1e-05;
const int matrix_tag = 123;
const int x_tag = 124;
const int b_tag = 125;

void create_matrix(double *matrix, int size, double mainData, double otherData) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i == j)
                *(matrix + i * n + j) = mainData;
            else
                *(matrix + i * n + j) = otherData;
        }
    }
}

void print_matrix(const double *matrix, int size1, int size2) {
    for (int i = 0; i < size1; ++i) {
        for (int j = 0; j < size2; ++j) {
            printf("%lf ", *(matrix + i * n + j));
        }
        printf("\n");
    }
}

double square_sum(const double *vector, int size) {
    double sum = 0.0;
    for (int i = 0; i < size; ++i) {
        sum += *(vector + i) * *(vector + i);
    }
    return sum;
}

void send_matrix_to_processes(double *matrix, int matrix_size, int process_count, int tag) {
    for (int i = 1; i < process_count; ++i) {
        if (i < n % process_count) {
            MPI_Send(matrix + i * n * (n / process_count + 1), n * (n / process_count + 1), MPI_DOUBLE, i, tag,
                     MPI_COMM_WORLD);
        } else {
            MPI_Send(matrix + (n % process_count) * n * (n / process_count + 1) + (i - 1) * n * (n / process_count),
                     n * (n / process_count),
                     MPI_DOUBLE, i, 123, MPI_COMM_WORLD);
        }
    }
}

void send_vector_to_processes(double *vector, int vector_size, int process_count, int tag) {
    for (int i = 1; i < process_count; ++i) {
        if (i < n % process_count) {
            MPI_Send(vector + i * (n / process_count + 1), (n / process_count + 1), MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
        } else {
            MPI_Send(vector + (n % process_count) * (n / process_count + 1) + (i - 1) * (n / process_count),
                     (n / process_count),
                     MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
        }
    }
}

void printVector(const double *vector, int size) {
    for (int i = 0; i < size; ++i) {
        printf("%lf ", *(vector + i));
    }
    printf("\n");
}

int main(int argc, char **argv) {
    int rank, size;
    double *matrix = NULL, *x = NULL, *b = NULL;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int rowsCount = (rank < n % size ? (n / size + 1) : (n / size));
    if (rank == 0) {
        matrix = malloc(sizeof(double) * n * n);
        x = calloc(n, sizeof(double));
        b = malloc(sizeof(double) * n);
        for (int i = 0; i < n + 1; ++i) {
            *(b + i) = n + 1;
        }
        create_matrix(matrix, n, 2.0, 1.0);
        send_matrix_to_processes(matrix, n, size, matrix_tag);
        send_vector_to_processes(x, n, size, x_tag);
        send_vector_to_processes(b, n, size, b_tag);
    } else {
        matrix = malloc(sizeof(double) * n * rowsCount);
        x = malloc(sizeof(double) * rowsCount);
        b = malloc(sizeof(double) * rowsCount);
        MPI_Recv(matrix, n * rowsCount, MPI_DOUBLE, 0, matrix_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(x, rowsCount, MPI_DOUBLE, 0, x_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b, rowsCount, MPI_DOUBLE, 0, b_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //print_matrix(matrix, rowsCount, n);
        printVector(b, rowsCount);

    }
    MPI_Finalize();
    return 0;
}
