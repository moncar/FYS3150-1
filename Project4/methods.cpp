#include "methods.h" //headerfile
#include "triDiagonalSolver.h" //header tridiag
#include <cmath> //acos(-1)=pi

//this can(and should...) be more generalized and compressed...

//call to constructor
Methods::Methods(double mtime, double mdt, double mdx) {
    time = mtime; //total time
    dt = mdt; //time-step
    dx = mdx; //step-length
    n = 1./dx + 1; //vector size
    alpha = dt/(dx*dx); //constant alpha as defined
    if(alpha < 0.5) {
        exit(EXIT_FAILURE);
    }

    pi = std::acos(-1); //acos(-1)=pi
    x = arma::linspace(0,1,n);
    uPrev.set_size(n); //vector for previous values
    uNew.set_size(n); //vector for new values
    
    a.set_size(n); a.fill((-1.)*alpha); //lower- and upper diagonal
    b1.set_size(n); b1.fill(1+2*alpha); //mid diagonal in implicit
    b2.set_size(n); b2.fill(2+2*alpha); //mid diagonal in crank-nicolson
}

void Methods::setInitial() {
    /* function for intializing */ 
    uPrev = 1-x; //initial value (steady-state)
    for(int i=1; i<100 ;i++) {
        /* set initial value */
        uPrev += (-2./(pi*i))*arma::sin(i*pi*x);
    }
    return;
}

void Methods::explicitEuler() {
    /* Method Explicit Forward Euler */

    uNew.fill(0); //empty vector for solution
    
    Methods::setInitial(); //initialize

    for(int t=1; t<=(time/dt) ;t++) {
        /* Loop though time */
        for(int i=1; i<n-1 ;i++){
            /* algorithm as defined */
            uNew(i) = alpha*uPrev(i+1) + (1-2*alpha)*uPrev(i) + alpha*uPrev(i-1);
        }
        uNew(0) = 1; uNew(n-1) = 0;
        uPrev = uNew; //update values
    }

    std::cout << uNew << std::endl; //visualize
    
    return;
}

void Methods::implicitEuler() {
    /* Method Implicit Backward Euler */ 
    
    TriDiagonalSolver TGS = TriDiagonalSolver(); //object TriDiagonalSolver

    Methods::setInitial(); //initialize

    arma::mat A(n,n, arma::fill::zeros);
    arma::vec ai(n-1); ai.fill(-alpha);
    arma::vec bi(n); bi.fill(1+2*alpha);

    A.diag(-1) = ai; A.diag(0) = bi; A.diag(1) = ai;

    /* the algorithm */
    for(int t=1; t<=(time/dt) ;t++) {
        /* Loop through time */
        uNew = TGS.solver(n,a,b1,a,uPrev); //solver in TGS
        //uNew = arma::solve(A,uPrev);
        uNew(0) = 1.; uNew(n-1) = 0.;
        uPrev = uNew; //update values
    }

    std::cout << uNew << std::endl; //visualize
    
    return;
}

void Methods::implicitCrankNicolson() {
    /* Method Implicit Crank-Nicolson */
    
    arma::vec uTemp(n, arma::fill::zeros); //empty vector for temporary values
    
    TriDiagonalSolver TGS = TriDiagonalSolver(); //object TriDiagonalSolver

    Methods::setInitial(); //initialize

    //std::cout << uPrev << std::endl;
    for(int t=0; t<=(time/dt) ;t++) {
        /* Loop through time */
        for (int i=1; i<n-1 ;i++) {
            /* temporary with Explicit */
            uTemp(i) = alpha*uPrev(i-1) + (2-2*alpha)*uPrev(i) + alpha*uPrev(i+1);
        }
        uNew = TGS.solver(n,a,b2,a,uTemp); //solver in TGS
        uNew(0) = 1; uNew(n-1) = 0; //set boundary
        uPrev = uNew; //update values
    }
    
    std::cout << uNew << std::endl; //visualize

    return;
} //end function implicitCrankNicolson
