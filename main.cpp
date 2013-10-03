using namespace std;

#include <iostream>
#include <armadillo>
#include <stdio.h>
#include <time.h>
#include "lib.h"
using namespace arma;

FILE * f;

int main(int argc, char** argv)
  {
    int N = 200;
    int i, j;
    double rho_min = 0.0, rho_max = 5.0;
    double omega = 1; // 0.5; 1.0; 5.0;
    double **A, ** Z, *d, *rho, *e, h = (rho_max-rho_min)/(N+1);
    clock_t t;
    d = new double[N];
    e = new double[N];
    rho = new double[N];


    for(i = 0; i < N; i++) 
      rho[i] = (i+1)*h;

    A = (double **) matrix(N, N, sizeof(double));
    Z = (double **) matrix(N, N, sizeof(double));

    // Initialize tridiagonal array
    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
            if(i == j)
            {
                A[i][j] = 2.0/(h*h) + omega*omega*rho[i]*rho[i] + 1/rho[i];
                d[i] = A[i][j];
                Z[i][j] = 1.0;
            }
            else if((j-1 == i) || (i-1 == j) )
                {
                   A[i][j] = -1.0/(h*h);
                   e[i] = A[i][j];

                }
                else
            {
                A[i][j] = 0.0; Z[i][j] = 0.0;
            }

    tqli(d, e, N, Z);
    f = fopen("log.txt", "w");
    for(i = 1; i < N; i++)
        fprintf(f," %f %f\n", rho[i], Z[i][(int)round(d[i])]);

    delete [] d;
    delete [] e;
    free_matrix((void **) A);
    free_matrix((void **) Z);
    fclose(f);
    return 0;

}
