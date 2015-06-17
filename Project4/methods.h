#ifndef METHODS_H
#define METHODS_H

#include <armadillo>

class Methods {
    public:
        int n;
        double time;
        double dt;
        double dx;
        double alpha;
        double pi;
        
        arma::vec x;
        arma::vec uPrev;
        arma::vec uNew;

        arma::vec a;
        arma::vec b1;
        arma::vec b2;
        
        Methods(double, double, double);

        void setInitial();
        void explicitEuler();
        void implicitEuler();
        void implicitCrankNicolson();
};

#endif //METHODS_H
