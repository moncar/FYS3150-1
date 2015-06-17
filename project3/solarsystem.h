#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include "body.h"
#include <vector>

class Solarsystem {
public:
    std::vector<Body> objects;
    double kinetic_energy;
    double potential_energy;
    double total_energy();
    arma::vec angular_momentum;

    Solarsystem();
    void addBody(Body newObject);
    void calculate_force_and_energy();
    int number_of_bodies();
};

#endif // SOLARSYSTEM_H

