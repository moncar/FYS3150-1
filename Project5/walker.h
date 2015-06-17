#ifndef WALKER_H
#define WALKER_H

#include <vector> //vector
#include <string> //string
#include <armadillo> //arma mat
#include "randomNumbers.h" //header for randomNumber

class Walker {
    public:
        double walkProbability;
        int numberWalkers;
        int maxTrials;

        double dt;
        double l_0;
        double minDistance;
        double maxDistance;
        int positionSize;

        bool useGauss;

        RandomNumber *RN;

        Walker(int, double, double, double, double, double, RandomNumber*);

        void setParameters(double, double);
        void startGauss(bool);
        void mcSampling(std::vector<double>&, std::vector<double>&, std::vector<double>&);
        void mcSampling2(std::vector<double>&, std::vector<double>&, std::vector<double>&,
                std::vector<double>&, std::vector<double>&, std::vector<double>&,
                arma::mat&);
        void output1(const char*, const char*, 
                std::vector<double>, std::vector<double>, std::vector<double>);
        void output2(const char*, const char*, 
                std::vector<double>, std::vector<double>, std::vector<double>,
                std::vector<double>, std::vector<double>, std::vector<double>);
};

#endif //WALKER_H
