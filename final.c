#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

const int n = 10;
const int tag = 123;
const double t = 0.01;
const double eps = 1e-07;

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

void mpi_start(int *argc, char ***argv, int *rank, int *size) {
    MPI_Init(argc, argv);
    MPI_Comm_size(MPI_COMM_WORLD, size);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
}

void alloc_vectors(double **x, double **b, double **next_x, double **main_x, int **matrix_displs, int **rows_counts,
                   int **vector_displs, int **send_counts, int proc_size) {
    *x = calloc(n, sizeof(double));
    *b = malloc(n * sizeof(double));
    *next_x = malloc(n * sizeof(double));
    *matrix_displs = calloc(proc_size, sizeof(int));
    *rows_counts = malloc(proc_size * sizeof(int));
    *main_x = calloc(n, sizeof(double));
    *vector_displs = calloc(proc_size, sizeof(int));
    *send_counts = calloc(proc_size, sizeof(int));
}

void free_vectors(double **x, double **b, double **next_x, double **main_x, int **matrix_displs, int **rows_counts,
                  int **vector_displs, int **send_counts) {
    free(*x);
    free(*b);
    free(*next_x);
    free(*matrix_displs);
    free(*rows_counts);
    free(*main_x);
    free(*vector_displs);
    free(*send_counts);
}

void send_matrix_to_processes(double *matrix, int *matrix_displs, int *rows_counts, int proc_size) {
    for (int i = 1; i < proc_size; ++i) {
        MPI_Send(matrix + matrix_displs[i], n * rows_counts[i], MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
    }
}

void
init_vectors(double *b, int *matrix_displs, int *rows_counts, int *vector_displs, int *send_counts, int proc_size) {
    for (int i = 0; i < n; ++i) {
        *(b + i) = n + 1;
    }
    for (int i = 0; i < proc_size; ++i) {
        *(rows_counts + i) = (i < n % proc_size) ? (n / proc_size) + 1 : (n / proc_size);
        *(send_counts + i) = n;
    }
    for (int i = 0; i < proc_size; ++i) {
        *(matrix_displs + i) = *(matrix_displs + i - 1) + n * *(rows_counts + i - 1);
        *(vector_displs + i) = *(matrix_displs + i) / n;
    }
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

void iteration(double *x, double *b, double *matrix, double *result, int rows_count, int column_count) {
    double *ax = malloc(sizeof(double) * rows_count);
    double *axB = malloc(sizeof(double) * rows_count);
    mult_matrix_vector(matrix, x, ax, rows_count, column_count);
    sub_vect(ax, b, axB, rows_count);
    mult_vec_digit(axB, t, rows_count);
    sub_vect(x, axB, result, rows_count);
    free(ax);
    free(axB);
}

void printVector(const double *vector, int size) {
    for (int i = 0; i < size; ++i) {
        printf("%lf ", *(vector + i));
    }
    printf("\n");
}

double vector_square_sum(const double *vector, int rows_count) {
    double result = 0.0;
    for (int i = 0; i < rows_count; ++i) {
        result += *(vector + i) * *(vector + i);
    }
    return result;
}


int main(int argc, char **argv) {
    int rank, size;
    int *matrix_displs, *rows_counts, *vector_displs, *send_counts;
    double *matrix = NULL, *x, *main_x, *b, *next_x;
    double local_square_sum, global_square_sum;
    mpi_start(&argc, &argv, &rank, &size);
    alloc_vectors(&x, &b, &next_x, &main_x, &matrix_displs, &rows_counts, &vector_displs, &send_counts, size);
    init_vectors(b, matrix_displs, rows_counts, vector_displs, send_counts, size);
    if (rank == 0) {
        matrix = malloc(sizeof(double) * n * n);
        create_matrix(matrix, n, n, 2.0, 1.0);
        send_matrix_to_processes(matrix, matrix_displs, rows_counts, size);
    } else {
        matrix = malloc(sizeof(double) * n * rows_counts[rank]);
        MPI_Recv(matrix, n * rows_counts[rank], MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    double *ax = malloc(sizeof(double) * rows_counts[rank]);
    double *axB = malloc(sizeof(double) * rows_counts[rank]);
    const double start_time_s = MPI_Wtime();
    do {
        iteration(x, b, matrix, next_x, rows_counts[rank], n);
        MPI_Gatherv(next_x, rows_counts[rank], MPI_DOUBLE, main_x, rows_counts, vector_displs, MPI_DOUBLE, 0,
                    MPI_COMM_WORLD);
        int *scatter_displs = calloc(size, sizeof(int));
        MPI_Scatterv(main_x, send_counts, scatter_displs, MPI_DOUBLE, x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        mult_matrix_vector(matrix, x, ax, rows_counts[rank], n);
        sub_vect(ax, b, axB, rows_counts[rank]);
        local_square_sum = vector_square_sum(axB, rows_counts[rank]);
        MPI_Reduce(&local_square_sum, &global_square_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Bcast(&global_square_sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        free(scatter_displs);

    } while (global_square_sum / vector_square_sum(b, n) >= eps * eps);
    const double end_time_s = MPI_Wtime();
    //printf("Process: %d , %f\n", rank, end_time_s - start_time_s);
    printVector(x, n);
    free_vectors(&x, &b, &next_x, &main_x, &matrix_displs, &rows_counts, &vector_displs, &send_counts);
    MPI_Finalize();
    return 0;
}