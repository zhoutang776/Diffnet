#include <R.h>
#include <cassert>
#include <iostream>
#define ARMA_DONT_USE_WRAPPER
#include "./include/armadillo"

using namespace std;
using namespace arma;

extern "C" {
mat soft(mat X, double lambda){
    mat Y = X;
    Y = (Y>lambda)%(Y-lambda)+(Y<(-lambda))%(Y+lambda);
    return Y;
}


inline double loss_func(const mat &Sigma_X, const mat &Sigma_Y, const mat &diff_Sigma, const mat &Delta){
    return 0.25*sum(sum(Delta % (Sigma_X*Delta*Sigma_Y))) + 0.25*sum(sum(Delta % (Sigma_Y*Delta*Sigma_X)))  - sum(sum(Delta % diff_Sigma));
}

inline double penalty_func(const mat &Delta, double lambda){
    return  lambda * sum(sum(abs(Delta)));
}


void symmetric_admm(double * pDelta0, double * pDelta3, double * pLambda0, double * pLambda3,
		 double * pSigmaX, double * pSigmaY, int * rho,
         double * pC1, double * pC2, double * pUx, double * pUy, double * tol, int * p, double * lambda, int * max_iter, int * iter)
{
    const mat SigmaX(pSigmaX, p[0], p[0], false);
    const mat SigmaY(pSigmaY, p[0], p[0], false);
    const mat diff_Sigma = SigmaX - SigmaY;

    const mat C1(pC1, p[0], p[0], false);
    const mat C2(pC2, p[0], p[0], false);

    const mat Ux(pUx, p[0], p[0], false);
    const mat Uy(pUy, p[0], p[0], false);


    mat A(p[0], p[0]), B(p[0], p[0]), C(p[0], p[0]);

    mat Delta1(pDelta0, p[0], p[0]), Delta2(pDelta0, p[0], p[0]);
    mat Delta3(pDelta3, p[0], p[0], false, true);
    mat temp(p[0], p[0]);

    mat Lambda1(pLambda0, p[0], p[0]), Lambda2(pLambda0, p[0], p[0]);
    mat Lambda3(pLambda3, p[0], p[0], false, true);


    
	iter[0] = 0;
    double f_old = 0, f = 0, err = 0;

    f = loss_func(SigmaX, SigmaY, diff_Sigma, Delta3) + penalty_func(Delta3, lambda[0]);

    while(iter[0] < max_iter[0]){

        A = SigmaX - SigmaY + rho[0]*Delta2 + rho[0]*Delta3 + Lambda3 - Lambda1;

        Delta1 = Uy *(C1 % (Uy.t()*A*Ux))* Ux.t();

        B = SigmaX - SigmaY + rho[0]*Delta1 + rho[0]*Delta3 + Lambda1 - Lambda2;

        Delta2 = Ux *(C2 % (Ux.t()*B*Uy))* Uy.t();

        C = (Lambda2/rho[0] - Lambda3/rho[0] + Delta1 + Delta2) / 2;

        Delta3 = soft(C, (*lambda / *rho) / 2);

		/*
        for(i=0;i<p[0];i++)
            for(j=0;j<p[0];j++){
                if(Delta3.at(i,j) == 0) Delta3.at(j,i)=0;
            }
        */
        temp = Delta3;
        Delta3 = (temp + temp.t())/2;
        

        f_old = f;
        f = loss_func(SigmaX, SigmaY, diff_Sigma, Delta3) + penalty_func(Delta3, lambda[0]);

        err = (f_old - f)/(f_old+1);
        err = abs(err);

        if (err < tol[0] && iter[0] !=0 ){
            break;
        }

        Lambda1 = Lambda1 + rho[0] * (Delta1 - Delta2);
        Lambda2 = Lambda2 + rho[0] * (Delta2 - Delta3);
        Lambda3 = Lambda3 + rho[0] * (Delta3 - Delta1);

        iter[0] = iter[0] + 1;
    }
    //if(m == max_iter[0]) cout<<"symmetric ADMM converge fails!!"<<endl;
    
}
}
