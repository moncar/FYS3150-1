#ifndef TRIDIAGONALSOLVER_H
#define TRIDIAGONALSOLVER_H

#include <armadillo>

class TriDiagonalSolver {
    public:
        TriDiagonalSolver();

        arma::vec solver(int n, arma::vec a, arma::vec b, arma::vec c, arma::vec RHS);
};

#endif //TRIDIAGONALSOLVER_H

