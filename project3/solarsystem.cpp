#include "solarsystem.h" //header file
#include <cmath> //std::acos(double) //acos(-1)=pi

Solarsystem::Solarsystem() {
    /* construct system */
}

void Solarsystem::addBody(Body newObject) {
    /* function for adding an object to system */
    objects.push_back(newObject);
}

void Solarsystem::calculate_force_and_energy() {
    /* function for calculating force, energy and angular momentum*/
    
    //set values to zero
    kinetic_energy = 0;
    potential_energy = 0;
    angular_momentum.zeros();

    double pi = std::acos(-1); //pi

    for(int i=0; i<number_of_bodies() ;i++) {
        /* loop through all bodies */
        Body &body1 = objects[i]; //initialize first body
        for(int j=0; j<number_of_bodies() ;j++) {
            if(j != i) {
                /* don't want to calculate forces twice */
                Body &body2 = objects[j]; //initialize second body
                arma::vec delta_R = body1.position - body2.position; //vector distance between bodies
                double dr = arma::norm(delta_R); //length of vectori above
          
                //calulate force
                body2.force += (4*pi*pi*body1.mass*body2.mass/(dr*dr*dr))*delta_R;
                
                //calulate potential energy -GM/r
                potential_energy += ((-1)*(-4*pi*pi)*body1.mass)/dr;
            }
        }
        
        //calculate kinetic energy 1/2mv^2
        double velocity = arma::norm(body1.velocity);
        kinetic_energy += (1/2.)*body1.mass*velocity*velocity;

        //caluclate angular momentum r x p = r x mv
        angular_momentum = arma::cross(body1.position,body1.mass*body1.velocity);
    }
    return;
}

int Solarsystem::number_of_bodies() {
    /* function returning number of bodies */
    return objects.size();
}

double Solarsystem::total_energy() {
    /* function calculating total energy of system */
    return kinetic_energy + potential_energy;
}
