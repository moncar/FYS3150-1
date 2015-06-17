#include "methods.h" //header file
#include <iostream> //cout
#include <cstdio> //fopen, fclose, fprintf

//call to constructor
Methods::Methods() {
}

void Methods::rk4_integrate(Solarsystem &system, double total_time, double dt) {
    /* method RungeKutta4 */

    int N = total_time/dt + 1; //number of values
   
    //create files
    FILE *outputfiles[10] = {
        fopen("sun.txt", "w"), fopen("mercury.txt", "w"), fopen("venus.txt", "w"), fopen("earth.txt", "w"), fopen("mars.txt", "w"), 
        fopen("jupiter.txt", "w"), fopen("saturn.txt", "w"), fopen("uranus.txt", "w"), fopen("neptune.txt", "w"), fopen("pluto.txt", "w")};
    
    for(int k=0; k<=N; k++) {
        for(int i=0; i<system.number_of_bodies(); i++) {
            Body &object = system.objects[i]; //initialize object
      
            //write to files
            std::fprintf(outputfiles[i], "%15f %15f %15f\n", object.position[0], object.position[1], object.position[2]);
       
            // RK4 algorithm //
            
            //previous values of position and velocity
            arma::vec rp = object.position; arma::vec vp = object.velocity;

            //reset the force on all bodies
            for(int j=0; j<system.number_of_bodies(); j++) {
                system.objects[j].resetForce();
            }

            //calulate new force, position and velocity(Eulers)
            system.calculate_force_and_energy();
            arma::vec K1 = vp*dt; 
            arma::vec L1 = (object.force/object.mass)*dt;
            object.velocity = vp + L1/2.;
            object.position = rp + K1/2.;
            
            //reset force on all bodies
            for(int j=0; j<system.number_of_bodies(); j++) {
                system.objects[j].resetForce();
            }

            //caluclate new force, position and velocity(at mid slope)
            system.calculate_force_and_energy();
            arma::vec K2 = object.velocity*dt;
            arma::vec L2 = (object.force/object.mass)*dt;
            object.velocity = vp + L2/2.;
            object.position = rp + K2/2.;
            
            //reset force on all bodies
            for(int j=0; j<system.number_of_bodies(); j++) {
                system.objects[j].resetForce();
            }

            //calculate new force, position and velocity(at second mid)
            system.calculate_force_and_energy();
            arma::vec K3 = object.velocity*dt;
            arma::vec L3 = (object.force/object.mass)*dt;
            object.velocity = vp + L3;
            object.position = rp + K3;

            //reset force on all bodies
            for(int j=0; j<system.number_of_bodies(); j++) {
                system.objects[j].resetForce();
            }

            //calculate final force, position and velocity(with previous values)
            system.calculate_force_and_energy();
            arma::vec K4 = object.velocity*dt;
            arma::vec L4 = (object.force/object.mass)*dt;

            //final position and velocity
            object.velocity = vp + (1./6)*(L1 + 2*L2 + 2*L3 + L4);
            object.position = rp + (1./6)*(K1 + 2*K2 + 2*K3 + K4);
        }
    }
    //close opened files after integration
    for(int j=0; j<sizeof outputfiles/sizeof (FILE*); j++) {
        fclose(outputfiles[j]);
    }
}

void Methods::verlet_integrate(Solarsystem &system, double total_time, double dt) {
    /* method velocity verlet */
    
    int N = total_time/dt + 1; //number of iterations (nr.o values)
    
    //create files 
    FILE *outputfiles[10] = {
        fopen("sun.txt", "w"), fopen("mercury.txt", "w"), fopen("venus", "w"), fopen("earth.txt", "w"), fopen("mars.txt", "w"), 
        fopen("jupiter.txt", "w"), fopen("saturn.txt", "w"), fopen("uranus", "w"), fopen("neptune.txt", "w"), fopen("pluto.txt", "w")};

    for(int k=0; k<=N; k++) {
        for(int i=0; i<system.number_of_bodies(); i++) {
            Body &object = system.objects[i]; //initialize body

            //write to file
            std::fprintf(outputfiles[i], "%15f %15f %15f\n", object.position[0], object.position[1], object.position[2]);
           
            //reset force on all bodies
            for(int j=0; j<system.number_of_bodies(); j++) {
                system.objects[j].resetForce();
            }

            // the algorithm //
            
            //calculate new force and predicted position
            system.calculate_force_and_energy();
            arma::vec K1 = 0.5*(object.force/object.mass)*dt;
            object.position += (object.velocity + K1)*dt;
            
            //reset force
            for(int j=0; j<system.number_of_bodies(); j++) {
                system.objects[j].resetForce();
            }

            //calulate new force with predicted value
            system.calculate_force_and_energy();
            arma::vec K2 = 0.5*(object.force/object.mass)*dt;
            
            //correct obtained value
            object.velocity += K1 + K2;
        }
    }
    //close opened files
    for(int j=0; j<sizeof outputfiles/sizeof (FILE*); j++) {
        fclose(outputfiles[j]);
    }
}
