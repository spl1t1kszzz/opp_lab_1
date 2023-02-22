#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

const int n = 10000;
const double t = 0.00001;
const double epsilon = 1e-05;
const int matrix_tag = 123;
const int x_tag = 124;
const int b_tag = 125;

void create_matrix(double *matrix, int row_count, int column_count, double mainData, double otherData) {
    for (int i = 0; i < row_count; ++i) {
        for (int j = 0; j < column_count; ++j) {
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

void send_matrix_to_processes(double *matrix, int process_count, int tag) {
    for (int i = 1; i < process_count; ++i) {
        if (i < n % process_count) {
            MPI_Send(matrix + i * n * (n / process_count + 1), n * (n / process_count + 1), MPI_DOUBLE, i, tag,
                     MPI_COMM_WORLD);
        } else {
            int j = (n % process_count == 0) ? i : (i - 1);
            MPI_Send(matrix + (n % process_count) * n * (n / process_count + 1) + j * n * (n / process_count),
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


void sub_vect(const double *left, const double *right, double *result, int size) {
    for (int i = 0; i < size; ++i) {
        *(result + i) = *(left + i) - *(right + i);
    }
}

void mult_vec_digit(double *vector, double digit, int size) {
    for (int i = 0; i < size; ++i) {
        *(vector + i) = *(vector + i) * digit;
    }
}

double vector_square_sum(const double *vector, int rows_count) {
    double result = 0.0;
    for (int i = 0; i < rows_count; ++i) {
        result += *(vector + i) * *(vector + i);
    }
    return result;
}


void mult_matrix_vector(const double *matrix, const double *vector, double *result, int rows_count, int column_count) {
    for (int i = 0; i < rows_count; ++i) {
        double part_sum = 0.0;
        for (int j = 0; j < column_count; ++j) {
            part_sum += *(matrix + i * n + j) * *(vector + j);
        }
        *(result + i) = part_sum;
    }
}

void iteration(double *x, double *b, double *matrix, double *result, double param, int rows_count, int column_count) {
    double *ax = malloc(sizeof(double) * rows_count);
    double *axB = malloc(sizeof(double) * rows_count);
    mult_matrix_vector(matrix, x, ax, rows_count, column_count);
    sub_vect(ax, b, axB, rows_count);
    mult_vec_digit(axB, t, rows_count);
    sub_vect(x, axB, result, rows_count);
    free(ax);
    free(axB);
}


void set_vector(const double *src, double *dest, int size) {
    for (int i = 0; i < size; ++i) {
        *(dest + i) = *(src + i);
    }
}

void init_vectors(double* b, int* r_counts_g, int* displs_g, int* counts_s, int size) {
    for (int i = 0; i < n; ++i) {
        *(b + i) = n + 1;

    }
    for (int i = 0; i < size; ++i) {
        *(r_counts_g + i) = (i < n % size ? (n / size + 1) : (n / size));
        *(counts_s + i) = n;
    }
    for (int i = 1; i < size; ++i) {
        *(displs_g + i) = *(displs_g + i - 1) + *(r_counts_g + i - 1);
    }
}


void alloc_vectors(double **x, int** r_counts_g, int** displs_g, double** b, double **next_x, double ** rbuf,
                   int** counts_s, int** displs_s, int size) {
    *x = calloc(n, sizeof(double));
    *next_x = calloc(n, sizeof(double));
    *rbuf = calloc(n, sizeof(double));
    *displs_g = calloc(size, sizeof(int));
    *r_counts_g = malloc(sizeof(int) * size);
    *b = malloc(sizeof(double) * n);
    *counts_s = malloc(sizeof(int) * size);
    *displs_s = calloc(size, sizeof(int));
}

void free_vectors(double **x, int** r_counts_g, int** displs_g, double** b, double **next_x, double ** rbuf,int** counts_s, int** displs_s) {
    free(*x);
    free(*rbuf);
    free(*r_counts_g);
    free(*displs_g);
    free(*b);
    free(*next_x);
    free(*counts_s);
    free(*displs_s);
}

void mpi_start(int* argc, char*** argv, int* rank, int* size) {
    MPI_Init(argc, argv);
    MPI_Comm_size(MPI_COMM_WORLD, size);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
}

int main(int argc, char **argv) {
    int rank, size;
    int *r_counts_g, *displs_g, *counts_s, *displs_s;
    double *matrix = NULL, *x, *b, *next_x, *rbuf;
    double b_norm = n * pow(n + 1, 2);
    mpi_start(&argc, &argv, &rank, &size);
    const double start_time_s = MPI_Wtime();
    alloc_vectors(&x, &r_counts_g, &displs_g, &b, &next_x, &rbuf,&counts_s, &displs_s, size);
    init_vectors(b, r_counts_g, displs_g, counts_s, size);
    int row_count = (rank < n % size ? (n / size + 1) : (n / size));
    if (rank == 0) {
        matrix = malloc(sizeof(double) * n * n);
        create_matrix(matrix, n, n, 2.0, 1.0);
        send_matrix_to_processes(matrix,size, matrix_tag);

    } else {
        matrix = malloc(sizeof(double) * n * row_count);
        MPI_Recv(matrix, n * row_count, MPI_DOUBLE, 0, matrix_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    double* ax = malloc(sizeof(double ) * row_count);
    double* axB = malloc(sizeof(double) * row_count);
    double local_square_sum = 0.0;
    double global_square_sum = 0.0;
    clock_t start = clock();
    do {
        iteration(x, b, matrix, next_x, t, row_count, n);
        MPI_Gatherv(next_x, row_count, MPI_DOUBLE, rbuf, r_counts_g, displs_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(rbuf, counts_s, displs_s, MPI_DOUBLE, x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        mult_matrix_vector(matrix, x, ax, row_count, n);
        sub_vect(ax, b, axB, row_count);
        local_square_sum = square_sum(axB, row_count);
        MPI_Reduce(&local_square_sum, &global_square_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Bcast(&global_square_sum,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    while (global_square_sum/ b_norm >= epsilon * epsilon);
    const double end_time_s = MPI_Wtime();
    printf("Process: %d , %f\n", rank, end_time_s - start_time_s);
    free_vectors(&x, &r_counts_g, &displs_g, &b, &next_x, &rbuf, &counts_s, &displs_s);
    MPI_Finalize();
    return 0;
}
