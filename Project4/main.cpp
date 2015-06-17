#include "methods.h" //header file
#include <string> //string (converter)
#include <cstdlib> //atof

int main(int argc, char* argv[]) {
    double time = std::atof(argv[1]); //time in seconds

    double dx = 1./100; //delta x
    double dt = (dx*dx)/2.; //delta t

    /* run functions */
    Methods meth = Methods(time, dt, dx); //object meth from class Methods
        
    std::string method = argv[2]; //method of choice

    //if-tests for which method to run
    if(method == "EE") meth.explicitEuler();
    if(method == "IE") meth.implicitEuler();
    if(method == "CN") meth.implicitCrankNicolson();

    return 0;
}
