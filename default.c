#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

const int n = 200;
const long double t = 0.0000001;
const long double epsilon = 1e-07;

void createMatrix(long double *matrix, int size, long double mainData, long double otherData) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i == j)
                *(matrix + i * n + j) = mainData;
            else
                *(matrix + i * n + j) = otherData;
        }
    }
}

long double norm(const long double *vector, int size) {
    long double norm = 0.0;
    for (int i = 0; i < size; ++i) {
        norm += *(vector + i) * *(vector + i);
    }
    return sqrt(norm);
}


void printMatrix(const long double *matrix, int size) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            printf("%Lf ", *(matrix + i * n + j));
        }
        printf("\n");
    }
}

void printVector(const long double *vector, int size) {
    for (int i = 0; i < size; ++i) {
        printf("%Lf ", *(vector + i));
    }
    printf("\n");
}



void mult_mat_vec(const long double *matrix, const long double *vector, long double *result, int size) {
    for (int i = 0; i < size; ++i) {
        long double partSum = 0;
        for (int j = 0; j < size; ++j) {
            partSum += *(matrix + i * n + j) * *(vector + j);
        }
        *(result + i) = partSum;
    }
}

void mult_vec_digit(long double *vector, long double digit, int size) {
    for (int i = 0; i < size; ++i) {
        *(vector + i) = *(vector + i) * digit;
    }
}

void sub_vect(const long double *left, const long double *right, long double *result, int size) {
    for (int i = 0; i < size; ++i) {
        *(result + i) = *(left + i) - *(right + i);
    }
}

void iteration(long double *x, long double *b, long double *matrix, long double param, int size, long double *result) {
    long double *ax = malloc(sizeof(long double) * size);
    long double *axB = malloc(sizeof(long double) * size);
    mult_mat_vec((const long double *) matrix, x, ax, size);
    sub_vect(ax, b, axB, size);
    mult_vec_digit(axB, param, size);
    sub_vect(x, axB, result, size);
    free(ax);
    free(axB);
}

bool crit(long double *matrix, long double *x, long double *b, long double param, int size) {
    long double *ax = malloc(sizeof(long double) * size);
    long double *axB = malloc(sizeof(long double) * size);
    long double axBNorm, bNorm;
    mult_mat_vec((const long double *) matrix, x, ax, size);
    sub_vect(ax, b, axB, size);
    axBNorm = norm(axB, size);
    bNorm = norm(b, size);
    free(ax);
    free(axB);
    return (axBNorm / bNorm) <= param;
}


void set_vector(const long double *src, long double *dest, int size) {
    for (int i = 0; i < size; ++i) {
        *(dest + i) = *(src + i);
    }
}



int main(int argc, char **argv) {

    long double *matrix = malloc(sizeof(long double) * n * n);
    createMatrix(matrix, n, 2.0, 1.0);
    long double *x = calloc(n, sizeof(long double));
    long double *b = malloc(sizeof(long double) * n);
    for (int i = 0; i < n; ++i) {
        *(b + i) = n + 1;
    }
    long double *nextX = calloc(n, sizeof(long double));
    set_vector(x, nextX, n);
    clock_t start = clock();
    while (!crit(matrix, nextX, b, epsilon, n)) {
        iteration(x, b, matrix, t, n, nextX);
        set_vector(nextX, x, n);

    }
    printf("%ld\n", (clock() - start) / CLOCKS_PER_SEC);
    //printVector(nextX,n);
    free(matrix);
    free(x);
    free(nextX);
    free(b);
    return 0;
}
