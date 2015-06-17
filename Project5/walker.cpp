#include "walker.h" //header
#include <cmath> //sqrt
#include <iostream> //cout, endl
#include <string> //string

//construct walker
Walker::Walker(int mnumberWalkers, double mwalkProbability, 
        double time, double mdt, double mminDistance, double mmaxDistance,
        RandomNumber *mRN) {
    minDistance = mminDistance;
    maxDistance = mmaxDistance;

    walkProbability = mwalkProbability;
    numberWalkers = mnumberWalkers;

    Walker::setParameters(time, mdt);

    RN = mRN;
} //end constructor Walker

void Walker::setParameters(double time, double mdt) {
    /* function for initializing parameters */
    dt = mdt; //set time-step
    maxTrials = int (time/dt + 0.5); //set number of iterations
    l_0 = std::sqrt(2*dt); //l_0=sqrt(2*D*dt), defined as sqrt(5*dt) to keep explicit method stable
    positionSize = maxDistance/l_0 + 1; //set size of results
    
    return;
} //end function setParameters

void Walker::startGauss(bool gaussTest) {
    /* function for starting gauss */
    useGauss = gaussTest;

    return;
} //end function startGauss

void Walker::mcSampling(std::vector<double> &mxWalkerPosC, std::vector<double> &mxWalkerPosC2,
        std::vector<double> &xProbabilityPosition) {
    /* function mc-sampling for 1D */
    long seed = -1; //seed for random
    std::vector<double> position (numberWalkers,0); //position vector
    
    //loop over walkers
    for(int trial=0; trial<=maxTrials ;trial++) {
        //loop steps
        for(int walker=0; walker<=numberWalkers ;walker++) {
            if(position[walker]<minDistance || position[walker]>maxDistance) {
                /* keep walker leashed(boundary) */
                continue;
            } //end leashed

            double randomDistNum = RN->ran0(&seed); //random uniformly number
            double randomDistNumGauss; //initialize gauss random number
            if(useGauss == 1) {
                /* use gauss */
                randomDistNumGauss = RN->gauss(&seed); //gauss random number
            } else {
               randomDistNumGauss = 1.; //standard number
            }
            
            if(randomDistNum<walkProbability) {
                position[walker] += l_0*randomDistNumGauss;
            } else {
                position[walker] -= l_0*randomDistNumGauss;
            } //end random
            int index = (position[walker]-minDistance)/l_0 + 0.5;
            if(index < 0 || index >= positionSize) { 
                /* ignore if outside */
                continue;
            }
            //update cumulative positions
            mxWalkerPosC[trial] += position[walker];
            mxWalkerPosC2[trial] += position[walker]*position[walker];
            xProbabilityPosition[index] += 1;
        } //end walker
    } //end trial

    return;
} //end function mcSampling

void Walker::mcSampling2(std::vector<double> &mxWalkerPosC, std::vector<double> &mxWalkerPosC2, 
        std::vector<double> &xProbabilityPosition, 
        std::vector<double> &myWalkerPosC, std::vector<double> &myWalkerPosC2, 
        std::vector<double> &yProbabilityPosition,
        arma::mat &xyProbabilityPosition) {
    /* function mc-sampling for 2D */
    long seed = -1; //seed for random
    int xindex; int yindex;
    double randomDistNum;
    std::vector<double> xposition (numberWalkers,0); //position vector x
    std::vector<double> yposition (numberWalkers,0); //position vector y
    
    //loop over walkers
    for(int trial=0; trial<maxTrials ;trial++) {
        //loop steps
        for(int walker=0; walker<numberWalkers ;walker++) {
            if(xposition[walker]<minDistance || xposition[walker]>maxDistance) {
                /* keep x-walker leashed(boundary) */
                continue;
            } //end leashed x
            if(yposition[walker]<minDistance || yposition[walker]>maxDistance) {
                /* keep y-walker leashed(boundary) */
                continue;
            } //end leashed y

            randomDistNum = 2*RN->ran0(&seed); //random uniformly number

            if(randomDistNum<=walkProbability) {
                xposition[walker] += l_0;
            } else if(randomDistNum<=2*walkProbability) {
                xposition[walker] -= l_0;
            } else if(randomDistNum<=3*walkProbability) {
                yposition[walker] += l_0;
            } else {
                yposition[walker] -= l_0;
            } //end random

            xindex = (xposition[walker]-minDistance)/l_0 + 0.5;
            if(xindex < 0 || xindex >= positionSize) { 
                /* ignore outsiders x */
                continue;
            }
            
            yindex = (yposition[walker]-minDistance)/l_0 + 0.5;
            if(yindex < 0 || yindex >= positionSize) { 
                /* ignore outsiders y */
                continue;
            }
            
            //update cumulative positions x
            mxWalkerPosC[trial] += xposition[walker];
            mxWalkerPosC2[trial] += xposition[walker]*xposition[walker];
            xProbabilityPosition[xindex] += 1;

            //update cumulative positions y
            myWalkerPosC[trial] += yposition[walker];
            myWalkerPosC2[trial] += yposition[walker]*yposition[walker];
            yProbabilityPosition[yindex] += 1;

            xyProbabilityPosition(xindex,yindex) += 1.; //update 2D
        } //end walker
    } //end trial

    return;
} //end function mcSampling2

void Walker::output1(const char *name1, const char *name2,
        std::vector<double> mxWalkerPosC, std::vector<double> mxWalkerPosC2,
        std::vector<double> xProbabilityPosition) {
    /* function write data to file */
    FILE *outputfile1 = fopen(name1, "w");
    for(int i=0; i<=mxWalkerPosC.size() ;i++) {
        double xaverage = mxWalkerPosC[i]/maxTrials; 
        double xaverage2 = mxWalkerPosC2[i]/maxTrials;
        double xvariance = xaverage2 - xaverage*xaverage;
        std::fprintf(outputfile1, "%i %15f %15f\n", 
                i, xaverage, xvariance);
    }//end write1
    fclose(outputfile1); //close file1
    
    FILE *outputfile2 = fopen(name2, "w");
    for(int j=0; j<xProbabilityPosition.size() ;j++) {
        std::fprintf(outputfile2, "%15f %15f\n", 
                j*l_0, xProbabilityPosition[j]);
    } //end write2
    fclose(outputfile2); //close file2

    return;
} //end funcion output

void Walker::output2(const char *name1, const char *name2,
        std::vector<double> mxWalkerPosC, std::vector<double> mxWalkerPosC2,
        std::vector<double> xProbabilityPosition,
        std::vector<double> myWalkerPosC, std::vector<double> myWalkerPosC2,
        std::vector<double> yProbabilityPosition) {
    /* function write data to file */
    FILE *outputfile1 = fopen(name1, "w");
    for(int i=0; i<=mxWalkerPosC.size() ;i++) {
        double xaverage = mxWalkerPosC[i]/maxTrials; 
        double xaverage2 = mxWalkerPosC2[i]/maxTrials;
        double xvariance = xaverage2 - xaverage*xaverage;
        
        double yaverage = myWalkerPosC[i]/maxTrials; 
        double yaverage2 = myWalkerPosC2[i]/maxTrials;
        double yvariance = yaverage2 - yaverage*yaverage;
        std::fprintf(outputfile1, "%i %15f %15f %15f %15f\n", 
                i, xaverage, xvariance, yaverage, yvariance);
    }//end write1
    fclose(outputfile1); //close file1
    
    FILE *outputfile2 = fopen(name2, "w");
    for(int j=0; j<xProbabilityPosition.size() ;j++) {
        std::fprintf(outputfile2, "%15f %15f %15f\n", 
                j*l_0, xProbabilityPosition[j], yProbabilityPosition[j]);
    } //end write2
    fclose(outputfile2); //close file2

    return;
} //end funcion output
