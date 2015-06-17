#include "walker.h" //header walker (walker class)
#include "methods.h" //header file for explicit and implicit methods
#include "randomNumbers.h" //header file for ran0 and randomGauss
#include "cstdlib" //atof, atoi (convert command line arguments)
#include <iostream> //cout, endl (visualize)

int main(int argc, char *argv[]) {
    /* main function */
    double walkProbability = 0.5; //probability for direction (50/50)
    int numberWalkers = std::atoi(argv[1]); //number of particles
   
    double dt = std::atof(argv[2]); //time-step
    double minDistance = 0.; //lower boundary
    double maxDistance = 1.; //upper boundary

    double time = std::atof(argv[3]); //total time
  
    //initialize walker
    RandomNumber RN = RandomNumber();
    Walker wlk = Walker(numberWalkers, walkProbability, 
            time, dt, minDistance, maxDistance, &RN); //create walker

    Methods meth = Methods(&wlk); //create methods object (and send in walker)
   
    std::vector<double> xWalkerPosC (wlk.maxTrials,0); //position cumulative x
    std::vector<double> xWalkerPosC2 (wlk.maxTrials,0); //position squared cumulative x
    std::vector<double> xProbabilityPosition (wlk.positionSize,0); //position for walkers x
    
    std::vector<double> yWalkerPosC (wlk.maxTrials,0); //position cumulative y
    std::vector<double> yWalkerPosC2 (wlk.maxTrials,0); //position squared cumulative y
    std::vector<double> yProbabilityPosition (wlk.positionSize,0); //position for walkers y
  
    int gaussTest = std::atoi(argv[4]); //check for gauss through command
    bool sendTest; //initialize test variable
    if(gaussTest == 1) {
        /* check of command test is true */
        sendTest = true; //assign true to test
    } else {
        sendTest = false; //assign false to test
    }

    wlk.startGauss(sendTest); //use gauss (or not)
    
    std::string usage = argv[5];
    const char *filename1;
    const char *filename2;
    
    if(usage == "1D") {
        /* run 1D */
        wlk.mcSampling(xWalkerPosC, xWalkerPosC2, xProbabilityPosition); //run walk 1D
        const char *filename1 = "average_1D.txt";
        const char *filename2 = "test_1D.txt";

        wlk.output1(filename1, filename2, 
                xWalkerPosC, xWalkerPosC2, xProbabilityPosition); //run output 1D
    } else if(usage == "2D") {
        /* run 2D*/
        arma::mat xyProbabilityPosition(wlk.positionSize,wlk.positionSize, arma::fill::zeros);
        wlk.mcSampling2(xWalkerPosC, xWalkerPosC2, xProbabilityPosition,
                yWalkerPosC, yWalkerPosC2, yProbabilityPosition,
                xyProbabilityPosition); //run walk 2D
        
        const char *filename1 = "average_2D.txt";
        const char *filename2 = "test_2D.txt";
        const char *filename3 = "mc2D.txt";
       
        wlk.output2(filename1, filename2,
                xWalkerPosC, xWalkerPosC2, xProbabilityPosition, 
                yWalkerPosC, yWalkerPosC2, yProbabilityPosition); //run output 2D

        meth.output(filename3, xyProbabilityPosition); //output for 2D MC grid
    } else if(usage == "EE") {
        /* run explicit Euler */
        arma::mat solution = meth.explicitEuler2D(time);
        
        const char *filename1 = "expEuler2D.txt";
        meth.output(filename1, solution);
    } else if(usage == "IJ") {
        /* run implicit Jacobi */
        arma::mat solution = meth.implicitJacobi(time);

        const char *filename1 = "impJacobi2D.txt";
        meth.output(filename1, solution);
    } else {
        /* standard output for specification */
        std::cout << "Specify option, [options: 1D, 2D, EE, IJ]" << std::endl;
    } //end usage

    return 0;
} //end function main
