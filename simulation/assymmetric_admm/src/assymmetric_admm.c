#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "R.h"
#include <R.h>
#include <R_ext/BLAS.h>
void matrixTransport(double *a, double *b, int n){
    int i,j;
    for(i = 0;i < n;i++){
        for(j = 0; j < n ; j++){
            b[i*n+j] = a[j*n+i];
        }
    }
}

void matvecmult(double *a, double *b, int *nrow, int *ncol, double *result)
{
    double zero = 0.0;
    double one = 1.0;
    int ione = 1;
    F77_CALL(dgemv)("n", nrow, ncol, &one, a, nrow, b, &ione, &zero, result, &ione);
}

void matmatmult(double *a, double *b, int *nrowa, int *ncola, int *ncolb, double *c)
{
    double one = 1.0;
    double zero = 0.0;
    F77_CALL(dgemm)("n", "n", nrowa, ncolb, ncola, &one, a, nrowa, b, ncola, &zero, c, nrowa);
}



void assymmetric_admm(double * Delta0, double * Delta2, double * Lambda0, double * SigmaX, double * SigmaY, int * rho,
         double * C1,  double * Ux, double * Uy, double * tol, int * p, double * lambda)
{
	int dim;
	dim = *p;
	dim = pow(dim,2);

    double *A = (double*) malloc(dim*sizeof(double));
    double *B = (double*) malloc(dim*sizeof(double));
    double *C = (double*) malloc(dim*sizeof(double));
    double *Delta1 = (double*) malloc(dim*sizeof(double));
    double *Lambda1 = (double*) malloc(dim*sizeof(double));


    int i, j, m = 0;

    double *z1 = (double*) malloc(dim*sizeof(double));
	double *z3 = (double*) malloc(dim*sizeof(double));
    double *tz3 = (double*) malloc(dim*sizeof(double));
    double *A1 = (double*) malloc(dim*sizeof(double));
	double *A2 = (double*) malloc(dim*sizeof(double));
	double *A3 = (double*) malloc(dim*sizeof(double));
	double *A4 = (double*) malloc(dim*sizeof(double));
    double *tUy = (double*) malloc(dim*sizeof(double));
    double *tUx = (double*) malloc(dim*sizeof(double));

    double temp,a,b,c,d,dis1;
    dim=*p;
    // initialization
    for(i=0;i<dim;i++)
        for(j=0;j<dim;j++){
            Delta1[j*dim+i] = Delta0[j*dim+i];
            Delta2[j*dim+i] = Delta0[j*dim+i];
            Lambda1[j*dim+i] = Lambda0[j*dim+i];
        }
    matrixTransport(Uy, tUy, *p);
    matrixTransport(Ux, tUx, *p);

    while(m < 10000){
        // maybe it is C
        for(i=0; i<dim; i++)
            for(j=0;j<dim;j++){
                A[j*dim+i] = SigmaX[j*dim+i] - SigmaY[j*dim+i] + *rho * Delta2[j*dim+i] - Lambda1[j*dim+i];
            }
        for(i=0; i<dim; i++)
            for(j=0;j<dim;j++){
                A1[j*dim+i]=0;
                A2[j*dim+i]=0;
                A3[j*dim+i]=0;
                A4[j*dim+i]=0;
                z1[j*dim+i]=0;
            }
        //A2 = Ux' * A * Uy
        matmatmult(tUx, A, p, p, p, A1);
        matmatmult(A1, Uy, p, p, p, A2);
        // A3 = A2 .* C1
        for(i=0;i<dim;i++)
            for(j=0;j<dim;j++){
                A3[j*dim+i] = A2[j*dim+i] * C1[j*dim+i];
            }
        // z1 = Ux * A3 * Uy'
        matmatmult(Ux, A3, p, p, p, A4);
        matmatmult(A4, tUy, p, p, p, z1); // z1 is Delta1

        for(i=0;i<dim;i++)
            for(j=0;j<dim;j++){
                C[j*dim+i] = (Lambda1[j*dim+i] / *rho  + z1[j*dim+i] );
            }

        temp = (*lambda / *rho);

        for(i=0;i<dim;i++)
            for(j=0;j<dim;j++){
                z3[j*dim+i] = 0;
                if(C[j*dim+i] > temp)
                    z3[j*dim+i] = C[j*dim+i] - temp;
                if(C[j*dim+i] < (-temp))
                    z3[j*dim+i] = C[j*dim+i] + temp;
            } // z3 is Delta2
        
        

        for(i=0;i<dim;i++)
            for(j=0;j<dim;j++){
                if(z3[j*dim+i]==0){z3[j*dim+i]=0; z3[i*dim+j]=0;}
            }
		matrixTransport(z3, tz3, *p);
		
        for(i=0; i<dim; i++)
            for(j=0; j<dim; j++){
                z3[j*dim+i]=(z3[j*dim+i]+tz3[j*dim+i])/2;
            }
        // computing f-norm
        a=0;
        b=0;
        c=0;
        for(i=0; i<dim; i++)
            for(j=0;j<dim;j++){
                a += pow(Delta2[j*dim+i], 2);
                b += pow(z3[j*dim+i], 2);
                c += pow((Delta2[j*dim+i]-z3[j*dim+i]), 2);
            }

        a = sqrt(a);
        b = sqrt(b);
        c = sqrt(c);
        d = 1>a ? 1 : a;
        dis1 = c / (d>b ? d : b);

    
        if (dis1 < *tol){
            break;
        }
        for(i=0; i<dim; i++)
            for(j=0;j<dim;j++){
                Delta1[j*dim+i] = z1[j*dim+i];
                Delta2[j*dim+i] = z3[j*dim+i];                
                Lambda1[j*dim+i] = Lambda1[j*dim+i] + *rho * (Delta1[j*dim+i]-Delta2[j*dim+i]);

            }
        m = m + 1;
    }
    free(A);
    free(B);
    free(C);
    free(Delta1);
    free(Lambda1);

    free(z1);
	free(z3);
    free(tz3);
    free(A1);
    free(A2);
    free(A3);
  	free(A4);
    free(tUy);
    free(tUx);
}
