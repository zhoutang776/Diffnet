#include <R.h>
#include <cassert>
#include <iostream>
#define ARMA_DONT_USE_WRAPPER
#include "../include/armadillo"

#define BIG 1e20

using namespace std;
using namespace arma;
extern "C"
{
inline void shrinkage(mat& X, const mat& weighted_lambda) {
    X = sign(X) % clamp(abs(X) - weighted_lambda, 0, BIG);
}

inline void compute_weighted_lambda(const mat& Delta, mat& weighted_lambda, double lambda, double a){
    // assume Delta = abs(Delta)
    weighted_lambda = (Delta <= a*lambda) % (lambda - Delta/a);
}


void assymmetric_mcp(double *pSigma_X, double *pSigma_Y, int *p,
                 double *lambda, double *a, double *pX, double *pY, int *n_X, int *n_Y,
                 double *epsilon_X, double *epsilon_Y,
                 double *lip, double *stop_tol, int *max_iter, bool *verbose,
                 double *pDelta) {

    const mat Sigma_X(pSigma_X, p[0], p[0], false);
    const mat Sigma_Y(pSigma_Y, p[0], p[0], false);
    const mat diff_Sigma = Sigma_X - Sigma_Y;

    const mat X(pX, n_X[0], p[0], false); // must assume Sigma_X = trans(X) * X + epsilon_X[0]*eye(p[0], p[0]);
    const mat Y(pY, n_Y[0], p[0], false); // must assume Sigma_Y = trans(Y) * Y + epsilon_Y[0]*eye(p[0], p[0]);

    mat Delta(pDelta, p[0], p[0], false), Delta_inner_old(p[0], p[0]), Delta_outer_old(Delta), Delta_extra(p[0], p[0]);
    mat grad_L(p[0], p[0]);

    mat weighted_lambda = zeros<mat>(p[0], p[0]);

    double t, t_old;

    bool converge_flag = false;
    int outer_iter;

    double old_norm, new_norm, diff_norm;
    double inner_eta, outer_eta;
	
	// one step iteration
    for (outer_iter = 1; outer_iter <= 1; ++outer_iter){
        compute_weighted_lambda(abs(Delta_outer_old), weighted_lambda, lambda[0], a[0]);

        Delta_inner_old = Delta;
        // Delta_outer_old = Delta;
        t = 1; t_old = t;

        for (int inner_iter = 1; inner_iter <= max_iter[0]; ++inner_iter) {
            Delta_extra = Delta + (t_old - 1) / t * (Delta - Delta_inner_old);
			
	        if(n_X[0] > p[0] || n_Y[0] > p[0]){
            	grad_L = Sigma_X*Delta_extra*Sigma_Y - diff_Sigma;
		    }else{
		        if(epsilon_X[0] == 0 && epsilon_Y[0] == 0){
		            grad_L = trans(X) * (X * Delta_extra * trans(Y)) * Y - diff_Sigma;
		        }else{
		            grad_L = trans(X) * (X*Delta_extra*trans(Y)) * Y  +
		                     epsilon_Y[0] * trans(X) * (X*Delta_extra) +
		                     epsilon_X[0] * (Delta_extra*trans(Y)) * Y +
		                     epsilon_X[0] * epsilon_Y[0] * Delta_extra;
		            grad_L = grad_L  - diff_Sigma;
		        }
		    }		


            Delta_inner_old = Delta;
            Delta = Delta_extra - 1 / lip[0] * grad_L;
            shrinkage(Delta, weighted_lambda/lip[0]);

            t_old = t;
            t = (1 + sqrt(1 + 4 * pow(t, 2))) / 2;

            old_norm = norm(Delta_inner_old, "fro"); new_norm = norm(Delta, "fro");
            old_norm = 1 > old_norm ? 1 : old_norm;
            diff_norm = norm(Delta - Delta_inner_old, "fro");
            inner_eta = diff_norm / (old_norm > new_norm ? old_norm : new_norm);

            if (inner_eta < stop_tol[0]) {
            	converge_flag = 1;
                break;
            }
        }

        //old_norm = norm(Delta_outer_old, "fro"); new_norm = norm(Delta, "fro");
        //old_norm = 1 > old_norm ? 1 : old_norm;
        //diff_norm = norm(Delta - Delta_outer_old, "fro");
        //outer_eta = diff_norm / (old_norm > new_norm ? old_norm : new_norm);
        //if (outer_eta < outer_tol[0] && outer_iter > 1){
        //    break;
        //}

    }

    if (!converge_flag) {
        cout << "MCP not converges! lambda = " << lambda[0] << "  iter = " << max_iter[0]
             << "   eta =  " << inner_eta<< endl;
    }
    if(verbose[0] ){
        cout<<"lambda = "<<lambda[0]<<" iteration = "<<outer_iter<<endl;
    }

}

}
