mpicc default.c -o main -O3 -Wall -Wextra -Wpedantic --std=c99 -lm
mpecc -mpilog -O3 -o main main.c -Wall -Wpedantic --std=c99