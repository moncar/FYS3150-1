#include "randomNumbers.h"
#include <cmath>

RandomNumber::RandomNumber() {
}

double RandomNumber::ran0(long *idum) {
    /* function for uniform random number */
    const int a = 16807, m = 2147483647, q = 127773;
    const int r = 2836, MASK = 123459876;
    const double am = 1./m;
    long k; double ans;
   
    *idum ^= MASK;
    k = (*idum)/q;
    *idum = a*(*idum - k*q) - r*k;
    if(*idum < 0) {
        *idum += m;
    }

    ans = am*(*idum);
    *idum ^= MASK;

    return ans;
} //end function ran0

double RandomNumber::gauss(long *idum) {
    /* function for random number with gaussian distribution */
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;

    if(idum<0) {
        iset = 0;
    }
    if(iset==0) {
        do {
            v1 = 2.*RandomNumber::ran0(idum) - 1.;
            v2 = 2.*RandomNumber::ran0(idum) - 1.;
            rsq = v1*v1 + v2*v2;
        } while(rsq>=1. || rsq == 0.);
        fac = std::sqrt(-2.*std::log(rsq)/rsq);
        gset = v1*fac;
        iset = 1;
        return v2*fac;
    } else {
        iset = 0;
        return gset;
    }
} //end function gauss
