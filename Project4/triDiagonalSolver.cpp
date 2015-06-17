#include "triDiagonalSolver.h"

TriDiagonalSolver::TriDiagonalSolver() {
}

arma::vec TriDiagonalSolver::solver(int n, arma::vec a, arma::vec b, arma::vec c, arma::vec RHS) {
    arma::vec V(n, arma::fill::zeros); //solution
    arma::vec temp(n, arma::fill::zeros); //temporary
    double btemp; //temporary b-values

    //the algorithm
    btemp = b(1);
    V(1) = RHS(1)/btemp;
    for(int i=2; i<n ;i++) {
        temp(i) = c(i-1)/btemp;
        btemp = b(i) - a(i)*temp(i);
        V(i) = (RHS(i) - a(i)*V(i-1))/btemp;
    }
      
    for(int i=n-1; i>=1; i--) {
        V(i) -= temp(i+1)*V(i+1);
    }

    return V;
}
