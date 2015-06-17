#ifndef METHODS_H
#define METHODS_H

#include "solarsystem.h" //header file

class Methods {
    public:
        Methods();
        void rk4_integrate(Solarsystem &system, double total_time, double dt);
        void verlet_integrate(Solarsystem &system, double total_time, double dt);
};

#endif //METHODS_H
