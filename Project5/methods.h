#ifndef METHODS_H
#define METHODS_H

#include <armadillo>
#include "walker.h"

class Methods {
    public:
        Walker *wlk;
        int n;
        double h;
        double alpha;

        Methods(Walker*);

        arma::mat explicitEuler2D(double);
        arma::mat implicitJacobi(double time);
        void output(const char*, arma::mat);
};

#endif //METHODS_H
