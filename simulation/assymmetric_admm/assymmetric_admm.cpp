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
    return 0.5*sum(sum(Delta % (Sigma_X*Delta*Sigma_Y))) - sum(sum(Delta % diff_Sigma));
}

inline double penalty_func(const mat &Delta, double lambda){
    return  lambda * sum(sum(abs(Delta)));
}

void assymmetric_admm(double * pDelta, double * pLambda, double * pSigmaX, double * pSigmaY, double * rho,
                      double * pC,  double * pUx, double * pUy, double * tol, int * p, double * lambda, int *max_iter, int *iter)
{
    const mat Sigma_X(pSigmaX, p[0], p[0], false);
    const mat Sigma_Y(pSigmaY, p[0], p[0], false);
    const mat diff_Sigma = Sigma_X-Sigma_Y;

    const mat C(pC, p[0], p[0], false);
    const mat Ux(pUx, p[0], p[0], false);
    const mat Uy(pUy, p[0], p[0], false);

    mat Delta(pDelta, p[0], p[0], false, true);
    mat Lambda(pLambda, p[0], p[0], false, true);

    mat Z(Delta);


    iter[0] = 0;

    double err = 0, f_old = 0, f = 0;

    f = loss_func(Sigma_X, Sigma_Y, diff_Sigma, Z) + penalty_func(Z, lambda[0]);


    while((iter[0] < max_iter[0] && err > tol[0]) || iter[0] == 0){

        Delta = Ux * (C % (Ux.t() *(diff_Sigma + rho[0]*(Z - Lambda))* Uy) )* Uy.t();

        Z = soft(Delta + Lambda, lambda[0]/rho[0]);

        Lambda = Delta + Lambda - Z;

        //err = norm(Z - Z_old, "fro")/sqrt(p[0]);
        f_old = f;
        f = loss_func(Sigma_X, Sigma_Y, diff_Sigma, Z) + penalty_func(Z, lambda[0]);
        
        err = (f_old - f)/(f_old+1);
        err = std::fabs(err);

        iter[0] = iter[0] + 1;

    }
	//    cout<<err<<"   "<<tol[0]<<endl;

    //if(iter[0] == max_iter[0]) cout<<"assymmetric ADMM converge fails!"<<endl;

    Delta = Z;

}

};
