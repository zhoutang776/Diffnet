#include <R.h>
#include <cassert>
#include <iostream>
#define ARMA_DONT_USE_WRAPPER
#include "../include/armadillo"

#include <iomanip>

using namespace std;
using namespace arma;
extern "C"
{
mat soft(mat X, double lambda){
    mat Y=X;
    Y=(Y>lambda)%(Y-lambda)+(Y<(-lambda))%(Y+lambda);
    return Y;
}

void assymmetric_lasso(double *pSigma_X, double *pSigma_Y, int *p,
                       double *lambda, double *pX, double *pY, int *n_X, int *n_Y,
                       double *epsilon_X, double *epsilon_Y,
                       double *lip, double *stop_tol, int *max_iter, int *iter, bool *verbose,
                       double *pDelta) {

    const mat Sigma_X(pSigma_X, p[0], p[0], false);
    const mat Sigma_Y(pSigma_Y, p[0], p[0], false);
    const mat diff_Sigma = Sigma_X - Sigma_Y;

    const mat X(pX, n_X[0], p[0], false); // must assume Sigma_X = trans(X) * X + epsilon_X[0]*eye(p[0], p[0]);
    const mat Y(pY, n_Y[0], p[0], false); // must assume Sigma_Y = trans(Y) * Y + epsilon_Y[0]*eye(p[0], p[0]);

    mat Delta(pDelta, p[0], p[0], false), Delta_old = Delta, Delta_extra = Delta;


    double t = 1, t_old = t;


    double err = 0;
    iter[0] = 0;

    while(((err>stop_tol[0]) && (iter[0]<max_iter[0])) || iter[0] == 0){
        Delta_extra = Delta + (t_old - 1) / t * (Delta - Delta_old);

        Delta = soft(Delta_extra - (X.t()*(X*Delta_extra*Y.t())*Y-diff_Sigma)/lip[0], lambda[0]/lip[0]);

        t_old = t;
        t = (1 + sqrt(1 + 4 * pow(t, 2))) / 2;

        err = norm(Delta - Delta_old, "fro")/p[0];
        iter[0] += 1;
    }
    Delta_old = Delta;

    Delta = (Delta_old.t()+Delta_old) / 2;
    if (iter[0] == max_iter[0]) {
        cout << "LASSO not converges! lambda = " << lambda[0] << "  iter = " << max_iter[0]
             << "  error =  " << err << endl;
    }
    if(verbose[0]){
        cout<<"lambda = "<<lambda[0]<<" iteration = "<<iter[0]<<endl;
    }
}



}


//void aaa_assymmetric_lasso(double *pSigma_X, double *pSigma_Y, int *p,
//                       double *lambda, double *pX, double *pY, int *n_X, int *n_Y,
//                       double *epsilon_X, double *epsilon_Y,
//                       double *lip, double *stop_tol, int *max_iter, int *iter, bool *verbose,
//                       double *pDelta) {
//
//    const mat Sigma_X(pSigma_X, p[0], p[0], false);
//    const mat Sigma_Y(pSigma_Y, p[0], p[0], false);
//    const mat diff_Sigma = Sigma_X - Sigma_Y;
//
//    const mat X(pX, n_X[0], p[0], false); // must assume Sigma_X = trans(X) * X + epsilon_X[0]*eye(p[0], p[0]);
//    const mat Y(pY, n_Y[0], p[0], false); // must assume Sigma_Y = trans(Y) * Y + epsilon_Y[0]*eye(p[0], p[0]);
//
//    mat Delta(pDelta, p[0], p[0], false), Delta_old(Delta), Delta_extra(p[0], p[0]);
//    mat grad_L(p[0], p[0]);
//
//    //mat weighted_lambda = ones<mat>(p[0], p[0])*lambda[0];
//
//    double t = 1, t_old = t;
//
//    bool converge_flag = false;
//
//    double old_norm, new_norm, diff_norm, err;
//
//    for (iter[0] = 1; iter[0] <= max_iter[0]; ++iter[0]) {
//        Delta_extra = Delta + (t_old - 1) / t * (Delta - Delta_old);
//
//        if(n_X[0] >= p[0] || n_Y[0] >= p[0]){
//            grad_L = Sigma_X*Delta_extra*Sigma_Y - diff_Sigma;
//        }else{
//            if(epsilon_X[0] == 0 && epsilon_Y[0] == 0){
//                grad_L = trans(X) * (X * Delta_extra * trans(Y)) * Y - diff_Sigma;
//            }else{
//                grad_L = trans(X) * (X*Delta_extra*trans(Y)) * Y  +
//                         epsilon_Y[0] * trans(X) * (X*Delta_extra) +
//                         epsilon_X[0] * (Delta_extra*trans(Y)) * Y +
//                         epsilon_X[0] * epsilon_Y[0] * Delta_extra;
//                grad_L = grad_L  - diff_Sigma;
//            }
//        }
//
//
//        Delta_old = Delta;
//        Delta = soft(Delta_extra - 1 / lip[0] * grad_L, lambda[0]/lip[0]);
//
//        t_old = t;
//        t = (1 + sqrt(1 + 4 * pow(t, 2))) / 2;
//
////        old_norm = norm(Delta_old, "fro"); new_norm = norm(Delta, "fro");
////        diff_norm = norm(Delta - Delta_old, "fro");
////        old_norm = 1 > old_norm ? 1 : old_norm;
////        eta = diff_norm / (old_norm > new_norm ? old_norm : new_norm);
//        err = mean(mean(abs(Delta - Delta_old)));
//        if (err < stop_tol[0]) {
//            converge_flag = true;
//            break;
//        }
//    }
//    Delta = (Delta.t() + Delta)/2;
//    pDelta = Delta.memptr();
//
//    if (!converge_flag) {
//        cout << "LASSO not converges! lambda = " << lambda[0] << "  iter = " << max_iter[0]
//             << "  err =  " << err << endl;
//    }
//    if(verbose[0]){
//        cout<<"lambda = "<<lambda[0]<<" iteration = "<<iter[0]<<endl;
//    }
//}
