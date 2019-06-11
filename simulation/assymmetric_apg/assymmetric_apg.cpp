#include <R.h>
#include <cassert>
#include <iostream>
#define ARMA_DONT_USE_WRAPPER
#include "./include/armadillo"

#define BIG 1e20

using namespace std;
using namespace arma;
extern "C"
{
mat soft(mat X, mat Lambda){
    mat Y=X;
    Y=(Y>Lambda)%(Y-Lambda)+(Y<(-Lambda))%(Y+Lambda);
    return Y;
}

inline void compute_mcp_lambda(const mat& Delta, mat& weighted_lambda, double lambda, double a){
    // assume Delta = abs(Delta)
    weighted_lambda = (Delta <= a*lambda) % (lambda - Delta/a);
}

inline void compute_scad_lambda(const mat& Delta, mat& weighted_lambda, double lambda, double a){
    // assume Delta = abs(Delta)
    weighted_lambda = lambda * (
            (Delta<=lambda) + (Delta>lambda) % (clamp(a*lambda-Delta,0,BIG) / (a-1) / lambda)
    );
}

inline double small_p_loss_func(const mat &Sigma_X, const mat &Sigma_Y, const mat &diff_Sigma, const mat &Delta){
    return 0.5*sum(sum(Delta % (Sigma_X*Delta*Sigma_Y))) - sum(sum(Delta % diff_Sigma));
}

inline double big_p_loss_func(const mat &X, const mat &Y, const mat &diff_Sigma, const mat &Delta){
    // Delta must be symmetric!
    return 0.5*sum(sum(Delta % (X.t()*(X*Delta*Y.t())*Y) )) - sum(sum(Delta % diff_Sigma));
}

inline double penalty_func(const mat &Delta, const mat Lambda){
    return  sum(sum(abs(Lambda % Delta)));
}

void assymmetric_lasso(double *pSigma_X, double *pSigma_Y, int *p,
                  double *pLambda, double *pX, double *pY, int *n_X, int *n_Y,
                  double *epsilon_X, double *epsilon_Y,
                  double *lip, double *stop_tol, int *max_iter, int *iter, bool *verbose,
                  double *pDelta) {

    const mat Sigma_X(pSigma_X, p[0], p[0], false);
    const mat Sigma_Y(pSigma_Y, p[0], p[0], false);
	const mat diff_Sigma = Sigma_X - Sigma_Y;
    
	const mat Lambda(pLambda, p[0], p[0], false);
    
    const mat X(pX, n_X[0], p[0], false); // must assume Sigma_X = trans(X) * X + epsilon_X[0]*eye(p[0], p[0]);
    const mat Y(pY, n_Y[0], p[0], false); // must assume Sigma_Y = trans(Y) * Y + epsilon_Y[0]*eye(p[0], p[0]);

    mat Delta(pDelta, p[0], p[0], false, true), Delta_old(Delta), Delta_extra(Delta);
    mat grad_L(p[0], p[0]);

    
	iter[0] = 0;

    double err = 0, f_old = 0, f = 0;
	double t = 1, t_old = t;
    
	if(n_X[0] >= p[0] || n_Y[0] >= p[0]){
        f = small_p_loss_func(Sigma_X, Sigma_Y, diff_Sigma, Delta) + penalty_func(Delta, Lambda);
    }else{
        f = big_p_loss_func(X, Y, diff_Sigma, Delta) + penalty_func(Delta, Lambda);
    }
    

    while(((err>stop_tol[0]) && (iter[0]<max_iter[0])) || iter[0] == 0) {
        Delta_extra = Delta + (t_old - 1) / t * (Delta - Delta_old);
        Delta_old = Delta;

        // computing gradient
		// computing gradient
        if(n_X[0] >= p[0] || n_Y[0] >= p[0]){
            grad_L = Sigma_X*Delta_extra*Sigma_Y - diff_Sigma;
        }else{
            grad_L = trans(X)*(X*Delta_extra*trans(Y))*Y - diff_Sigma;
        }
        

        Delta = soft(Delta_extra - grad_L/lip[0], Lambda/lip[0]);

        f_old = f;
        if(n_X[0] >= p[0] || n_Y[0] >= p[0]){
            f = small_p_loss_func(Sigma_X, Sigma_Y, diff_Sigma, Delta) + penalty_func(Delta, Lambda);
        }else{
            f = big_p_loss_func(X, Y, diff_Sigma, Delta) + penalty_func(Delta, Lambda);
        }


        t_old = t;
        t = (1 + sqrt(1 + 4 * pow(t, 2))) / 2;

        // err = norm(Delta - Delta_old, "fro")/sqrt(p[0]);

        err = (f_old - f)/(f_old+1);
        err = std::fabs(err);

        iter[0] = iter[0] + 1;
    }
    
    //Delta_old = Delta;
    //Delta = (Delta_old + Delta_old.t())/2;

    if (iter[0] == max_iter[0])  {
        //cout << "assymmetric LASSO not converges! Iter = " << max_iter[0]<< "  error =  " << err << endl;
    }

}

void assymmetric_mcp(double *pSigma_X, double *pSigma_Y, int *p,
                 double *lambda, double *a, double *pX, double *pY, int *n_X, int *n_Y,
                 double *epsilon_X, double *epsilon_Y,
                 double *lip, double *stop_tol, int *max_iter, bool *verbose,
                 double *pDelta) {

    mat Delta(pDelta, p[0], p[0], false);
    mat weighted_lambda = ones<mat>(p[0], p[0])*lambda[0];
    int all_iter = 0, *iter = &all_iter;
    // one step iteration
    compute_mcp_lambda(abs(Delta), weighted_lambda, lambda[0], a[0]);

    assymmetric_lasso(pSigma_X, pSigma_Y, p, weighted_lambda.memptr(), pX, pY, n_X, n_Y,
                  epsilon_X, epsilon_Y,
                  lip, stop_tol, max_iter, iter, verbose,
                  pDelta);
    
    if(verbose[0] ){
        cout<<"lambda = "<<lambda[0]<<" iteration = "<<all_iter<<endl;
    }
}

void assymmetric_scad(double *pSigma_X, double *pSigma_Y, int *p,
                 double *lambda, double *a, double *pX, double *pY, int *n_X, int *n_Y,
                 double *epsilon_X, double *epsilon_Y,
                 double *lip, double *stop_tol, int *max_iter, bool *verbose,
                 double *pDelta) {

    mat Delta(pDelta, p[0], p[0], false);
    mat weighted_lambda = ones<mat>(p[0], p[0])*lambda[0];
    int all_iter = 0, *iter = &all_iter;
    // one step iteration
    compute_scad_lambda(abs(Delta), weighted_lambda, lambda[0], a[0]);

    assymmetric_lasso(pSigma_X, pSigma_Y, p, weighted_lambda.memptr(), pX, pY, n_X, n_Y,
                  epsilon_X, epsilon_Y,
                  lip, stop_tol, max_iter, iter, verbose,
                  pDelta);
    
    if(verbose[0] ){
        cout<<"lambda = "<<lambda[0]<<" iteration = "<<all_iter<<endl;
    }
}


}
