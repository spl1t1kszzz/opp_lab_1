mpicc final.c -o main -O3 -Wall -Wextra -Wpedantic --std=c99 -lm

mpecc -mpilog -O3 -o main final.c -Wall -Wpedantic --std=c99
