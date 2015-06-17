#ifndef BODY_H
#define BODY_H

#include<armadillo>

class Body {
    public:
        double mass;
        arma::vec position;
        arma::vec velocity;
        arma::vec force;
        Body(arma::vec pos, arma::vec vel, double mass_);
        void resetForce();
};

#endif //BODY_H
