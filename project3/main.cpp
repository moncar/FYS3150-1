#include "methods.h" //header file methods
#include <cmath> //acos(-1)=pi
#include <cstdlib> //atof
#include <string> //string

int main(int argc, char* argv[]) {
    double sun_mass     = 1.00000000e-00; //mass sun
    double mercury_mass = 1.65956463e-07; //mass mercury
    double venus_mass   = 2.44699613e-06; //mass venus
    double earth_mass   = 3.00245840e-06; //mass earth
    double mars_mass    = 3.00245840e-07; //mass mars
    double jupiter_mass = 9.54265748e-04; //mass jupiter
    double saturn_mass  = 2.85716656e-04; //mass saturn
    double uranus_mass  = 4.36430044e-05; //mass uranus
    double neptune_mass = 5.14855965e-05; //mass neptune
    double pluto_mass   = 6.58086572e-09; //mass pluto

    arma::vec mercury_position = {0.466697,0.,0.};
    arma::vec mercury_velocity = {0.,9.99077808,0.};

    arma::vec venus_position = {0.728213,0.,0.};
    arma::vec venus_velocity = {0.,7.38729464,0.};

    arma::vec earth_position = {1.,0,0.}; //initial position earth
    arma::vec earth_velocity = {0.,6.28194273,0.}; //initial velocity earth

    arma::vec mars_position {1.6660,0.,0.};
    arma::vec mars_velocity {0.,5.07892327,0.};

    arma::vec jupiter_position = {5.458104,0.,0.};
    arma::vec jupiter_velocity = {0.,2.75705142,0.};

    arma::vec saturn_position = {10.11595804,0.,0.};
    arma::vec saturn_velocity = {0.,2.04405725,0.};

    arma::vec uranus_position = {20.095371,0.,0.};
    arma::vec uranus_velocity = {0.,1.43442614,0.};

    arma::vec neptune_position = {30.331855,0.,0.};
    arma::vec neptune_velocity = {0.,1.14543146,0.};

    arma::vec pluto_position = {39.264,0.,0.};
    arma::vec pluto_velocity = {0.,0.991441599,0.};

    arma::vec velocities[9] = {mercury_velocity, venus_velocity, earth_velocity, mars_velocity, 
        jupiter_velocity, saturn_velocity, neptune_velocity, uranus_velocity, pluto_velocity};

    double masses[9] = {mercury_mass, venus_mass, earth_mass, mars_mass,
        jupiter_mass, saturn_mass, uranus_mass, neptune_mass, pluto_mass};

    double s_vel = 0;

    for(int i=0; i<9 ;i++) {
        s_vel += masses[i]*arma::norm(velocities[i]);
    }

    arma::vec sun_position = {0.,0.,0.}; //initial position sun
    arma::vec sun_velocity = {0.,s_vel,0.}; //initial position sun

    Solarsystem mySolarsystem; //initialize solarsystem

    Body sun    (sun_position, sun_velocity, sun_mass);             //create sun
    Body mercury(mercury_position, mercury_velocity, mercury_mass); //create mercury
    Body venus  (venus_position, venus_velocity, venus_mass);       //create venus
    Body earth  (earth_position, earth_velocity, earth_mass);       //create earth
    Body mars   (mars_position, mars_velocity, mars_mass);          //create mars
    Body jupiter(jupiter_position, jupiter_velocity, jupiter_mass); //create jupiter
    Body saturn (saturn_position, saturn_velocity, saturn_mass);    //create saturn
    Body uranus (uranus_position, uranus_velocity, uranus_mass);    //create uranus
    Body neptune(neptune_position, neptune_velocity, neptune_mass); //create neptune
    Body pluto  (pluto_position, pluto_velocity, pluto_mass);       //create pluto

    mySolarsystem.addBody(sun); //add sun to system
    mySolarsystem.addBody(mercury); //add mercury to system
    mySolarsystem.addBody(venus); //add venus to system
    mySolarsystem.addBody(earth); //add earth to system
    mySolarsystem.addBody(mars); //add mars to system
    mySolarsystem.addBody(jupiter); //add jupiter to system
    mySolarsystem.addBody(saturn); //add saturn to system
    mySolarsystem.addBody(uranus); //add uranus to system
    mySolarsystem.addBody(neptune); //add neptune to system
    mySolarsystem.addBody(pluto); //add pluto to system

    Methods methods = Methods(); //initialize methods

    double T  = std::atof(argv[1]); //total time
    double dt = std::atof(argv[2]); //time-step
    std::string method = argv[3]; //method choice

    if(method == "verlet") {
        /* run verlet */
        methods.verlet_integrate(mySolarsystem, T, dt); //call integrator with made solarsystem
    }
    else if(method == "rk4") {
        /* run rk4 */
        methods.rk4_integrate(mySolarsystem, T, dt);
    }
    else {
        /* standard */
        std::cout << "Specify method" << std::endl;
    }

    return 0;
}
