#include "body.h" //header file

Body::Body(arma::vec pos, arma::vec vel, double mass_) {
    /* construct body */
    position = pos; //position
    velocity = vel; //velocity
    mass = mass_; //mass
}

void Body::resetForce() {
    /* function for reseting force on body */
    force.zeros(3);
}
