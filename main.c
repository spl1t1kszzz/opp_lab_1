#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

const int n = 5;
const double t = 0.01;
const double epsilon = 1e-05;

void createMatrix(double **matrix, int size, double mainData, double otherData) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i == j) {
                *(*(matrix + i) + j) = mainData;
            } else {
                *(*(matrix + i) + j) = otherData;
            }
        }
    }
}

double norm(const double *vector, int size) {
    double norm = 0.0;
    for (int i = 0; i < size; ++i) {
        norm += *(vector + i) * *(vector + i);
    }
    return sqrt(norm);
}


void printMatrix(const double **matrix, int size) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            printf("%lf ", *(*(matrix + i) + j));
        }
        printf("\n");
    }
}

void printVector(const double *vector, int size) {
    for (int i = 0; i < size; ++i) {
        printf("%lf ", *(vector + i));
    }
    printf("\n");
}


double **allocMatrix() {
    double **matrix = (double **) malloc(sizeof(double *) * n);
    for (int i = 0; i < n; ++i) {
        *(matrix + i) = (double *) malloc(sizeof(double) * n);
    }
    return matrix;
}

void mult_mat_vec(const double **matrix, const double *vector, double *result, int size) {
    for (int i = 0; i < size; ++i) {
        double partSum = 0;
        for (int j = 0; j < size; ++j) {
            partSum += *(*(matrix + i) + j) * *(vector + j);
        }
        *(result + i) = partSum;
    }

}

void mult_vec_digit(double *vector, double digit, int size) {
    for (int i = 0; i < size; ++i) {
        *(vector + i) = *(vector + i) * digit;
    }
}

void sub_vect(const double *left, const double *right, double *result, int size) {
    for (int i = 0; i < size; ++i) {
        *(result + i) = *(left + i) - *(right);
    }
}

void iteration(double *vect, double *b, double **matrix, double param, int size, double *result) {
    double *ax = malloc(sizeof(double) * size);
    double *axB = malloc(sizeof(double) * size);
    mult_mat_vec((const double **) matrix, vect, ax, size);
    sub_vect(ax, b, axB, size);
    mult_vec_digit(axB, param, size);
    sub_vect(vect, axB, result, size);
    free(ax);
    free(axB);
}

bool crit(double **matrix, double *vect, double *b, double param, int size) {
    double *ax = malloc(sizeof(double) * size);
    double *axB = malloc(sizeof(double) * size);
    double axBNorm = 0, bNorm = 0;
    mult_mat_vec(matrix, vect, ax, size);
    sub_vect(ax, b, axB, size);
    axBNorm = norm(axB, size);
    bNorm = norm(b, size);
    free(ax);
    free(axB);
    return (axBNorm / bNorm) < param;
}

void freeMatrix(double **matrix) {
    for (int i = 0; i < n; ++i) {
        free(*(matrix + i));
    }
    free(matrix);
}

int main(/*int argc, char **argv*/) {
    double **matrix = allocMatrix();
    double *x = calloc(n, sizeof(double));
    double *b = malloc(sizeof(double) * n);
    for (int i = 0; i < n; ++i) {
        *(b + i) = n + 1;
    }
    double* nextX = calloc(n, sizeof(double));
    for (int i = 0; i < n; ++i) {
        *(nextX + i) = *(x + i);
    }
    createMatrix(matrix, n, 2.0, 1.0);
    while (!crit(matrix, nextX, b, epsilon, n)) {
        iteration(x, b, matrix, t, n, nextX);
        for (int i = 0; i < n; ++i) {
            *(x + i) = *(nextX + i);
        }
    }
    printVector(nextX, n);
    freeMatrix(matrix);
    free(x);
    free(nextX);
    free(b);
    return 0;
}
