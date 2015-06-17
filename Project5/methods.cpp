#include "methods.h" //header file
#include <cmath> //exp
#include <cstdlib> //exit

Methods::Methods(Walker *mwlk) {
    wlk = mwlk; //pull parameters from walker
    n = wlk->positionSize; //grid size
    h = wlk->maxDistance/(n+1); //step constant as defined
    alpha = wlk->dt/(h*h); //constant as defined
}

arma::mat Methods::explicitEuler2D(double time) {
    /* methods explicit euler (in 2D) */
    arma::mat uPrev(n,n, arma::fill::zeros); //empty vector for previous values
    arma::mat uNew(n,n, arma::fill::zeros); //empty vector for new values
    
    //set initial value
    for(int j=0; j<n ;j++) {
        for(int i=0; i<n ;i++) {
            uPrev(i,j) = (1-j*h)*std::exp(i*h); //steady-state
        } //end j
    } //end i
   
    
    for(int t=1; t<=(time/wlk->dt) ;t++) {
        /* calculate new values */
        for(int i=1; i<n-1 ;i++) {
            for(int j=1; j<n-1 ;j++) {
                uNew(i,j) = uPrev(i,j) + alpha*(uPrev(i+1,j)+uPrev(i-1,j)
                        +uPrev(i,j+1)+uPrev(i,j-1)-4*uPrev(i,j));
            } //end j
        } //end i
        
        for(int k=0; k<n ;k++) {
            /* set boundary conditions (assume x=y=[0,1]) */
            double kn   = float(k)/n;
            uNew(0,k)   = (1-kn)*std::exp(time);
            uNew(n-1,k) = (1-kn)*std::exp(1+time);
            uNew(k,0)   = std::exp(kn+time);
            uNew(k,n-1) = 0;
        } //end k

        uPrev = uNew; //update value
    } //end t
    
    return uNew;
} //end function explicitEuler2D

arma::mat Methods::implicitJacobi(double time) {
    /* method implicit Jacobi */
    arma::mat uPrev(n,n, arma::fill::zeros); //empty matrix for previous values
    arma::mat uNew(n,n, arma::fill::zeros); //empty matrix for new values
    arma::mat temp(n,n, arma::fill::zeros); //empty matrix for temporary values of uPrev
    
    double threshold = 0.000001; //convergance threshold
    double diff = threshold + 1; //difference of uPrev and temp

    //set initial value
    for(int j=0; j<n ;j++) {
        for(int i=0; i<n ;i++) {
            uPrev(i,j) = (1-j*h)*std::exp(i*h); //steady-state
        } //end j
    } //end i

    temp = uPrev; //set intial guess as previous value
   
    while(diff > threshold) {
        /* calculate new values */
        for(int i=1; i<n-1 ;i++) {
            for(int j=1; j<n-1 ;j++) {
                uNew(i,j) = (1./(1+4*alpha))*(uPrev(i,j) + 
                        alpha*(temp(i+1,j)+temp(i-1,j)+temp(i,j+1)+temp(i,j-1)));
            } //end j
        } //end i

        for(int k=0; k<n ;k++) {
            /* set boundary conditions (assume x=y=[0,1]) */
            double kn   = float(k)/n;
            uNew(0,k)   = (1-kn)*std::exp(time);
            uNew(n-1,k) = (1-kn)*std::exp(1+time);
            uNew(k,0)   = std::exp(kn+time);
            uNew(k,n-1) = 0;
        } //end k
        
        diff = 0; //reset diff
        for(int i=0; i<n ;i++) {
            /* calculate difference */
            for(int j=0; j<n ;j++) {
                diff += std::abs(temp(i,j) - uNew(i,j));
            } //end i
        } //end j
        diff /= n*n; //assuming our matrix is square
        
        temp = uPrev;
        uPrev = uNew; //keep temp
    } //end while

    return uNew;
} //end function implicitJacobi

void Methods::output(const char *filename, arma::mat results) {
    /* function for ouput to file */

    FILE *outputfile = fopen(filename, "w");
    std::fprintf(outputfile, "%i\n", n);
    for(int i=0; i<n ;i++) {
        for(int j=0; j<n ;j++) {
            std::fprintf(outputfile, "%15f\n", results(i,j));
        } //end i
    } //end j

    fclose(outputfile); //close file2
    return;
} //end function output
