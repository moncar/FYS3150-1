#include <iostream> //cout(show numbers)
#include <iomanip> //setprecision(decimal points)
#include <armadillo> //matrix and vector functions(mat, vec etc.)
#include <cstdlib> //atoi(convert input args), system(commandline input)
#include <cmath> //math functions(sqrt etc.)
#include <chrono> // timer(high_precision_clock)

using namespace arma;
using namespace std::chrono;


void print_and_sort(int number_of_values, vec eigenvalues) {
    /* function for sorting and printing eigenvalues */
    vec eigenvalues_sorted = sort(eigenvalues); //sort lowest-highest
    for(int i=0; i<number_of_values ;i++) {
        std::cout << eigenvalues_sorted[i] << std::endl;
    }
    return;
}

//Jacobi method as described(edited for armadillo matrices)
double maxoffdiag(int n, mat &A, int* k, int* l) {
    double max = 0.0;

    for(int i=0; i<n ;i++) {
        for(int j=i+1; j<n ;j++) {
            if(fabs(A(i,j)) > max) {
                max = fabs(A(i,j));
                *l = i; *k = j;
            }
        }
    }

    return max;
}

vec rotate(int n, mat &A, mat &R, int k, int l) {
    double s,c;
    if(A(k,l) != 0.0) {
        double t,tau;
        tau = (A(l,l) - A(k,k))/(2*A(k,l));
        if(tau > 0) {
            t = 1.0/(tau + std::sqrt(1. + tau*tau));
        } else {
            t = -1.0/(-tau + std::sqrt(1. + tau*tau));
        }
        c = 1./sqrt(1+t*t);
        s = c*t;
    } else {
        c = 1.0;
        s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    //changing matrix elements with indices k and l
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0; A(l,k) = 0.0;
    //change remaining elements
    for(int i=0; i<n; i++) {
        if(i != k && i != l) {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
            //std::cout <<
        }
        //compute new eigenvectors
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
    return A.diag();
}

void jacobi(int n, mat A, mat R) {
    int k,l;
    double tol = 1.0e-10;
    double max_it = n*n*n;
    int iterations = 0;
    double max_offdiag = maxoffdiag(n,A,&k,&l);
    
    vec eig;
    while(fabs(max_offdiag) >= tol && (double) iterations < max_it) {
        max_offdiag = maxoffdiag(n,A,&k,&l);
        eig = rotate(n,A,R,k,l);
        iterations++;
    }

    std::cout << std::setprecision(6) << std::fixed; //fix number of decimal points
    std::cout << "Eigenvalues-jacobi:" << std::endl;
    print_and_sort(5,eig);
    std::cout << "Iteratons:" << iterations << std::endl;
    
    return;
}
//end jacobi

void write_values_to_file(const char* filename, vec eigenvalues, 
        vec eigenvector1, vec eigenvector2, vec eigenvector3, 
        double rho_max, double rho_min, vec rho, int n) {
    /* create and write eigenvalues to seperate file.*/
    std::ofstream outputfile; //create filetype
    outputfile.open(filename); //create(and/or open) file

    outputfile << "rho_max =" << rho_max << '\n' 
        << "rho_min =" << rho_min << '\n' 
        << "n =" << n << '\n' 
        << eigenvalues 
        << eigenvector1 << eigenvector2 << eigenvector3 
        << rho;
    outputfile.close();
    
    return;
}

int main(int argc, char* argv[]) {
    if(argc < 2) std::cout << "USAGE: ./new_prog 'n' 'filname'" 
        << std::endl; //make sure n is included

    //constants defined in text
    int n_step = atoi(argv[1]); int n = n_step-1;
    double rho_max = 5.; double rho_min = 0.;
    double h = (rho_max - rho_min)/n_step;

    //values for rho and potensial V
    vec rho(n, fill::zeros); vec V(n, fill::zeros); 
    for(int i=0; i<n; i++) {
        rho(i) = rho_min + (i+1)*h;
        V(i) = rho(i)*rho(i);
    }

    //diagonal vector
    vec d(n); d.fill(2./(h*h));
    for(int i=0; i<n; i++) {
        d(i) += V(i);
    }    

    vec e(n-1); e.fill(-1./(h*h)); //Lower and Upper diagonal in A

    //create tridiagonal matrix A for eigenvalues
    mat A(n,n, fill::zeros);
    A.diag() = d; A.diag(-1) = e; A.diag(1) = e;
    
    //create matrix R for eigenvectors
    mat R(n,n, fill::zeros); R.diag().fill(1.);
    
    //run jacobi-method with timer
    auto start = high_resolution_clock::now();
    
    jacobi(n,A,R);

    auto finish = high_resolution_clock::now(); 
    std::cout << "time(jacobi):" 
        << duration_cast<seconds>(finish - start).count() << "s" 
        << std::endl;
    std::cout << "\n";
    
    //eigenvalue solver in armadillo with timer
    auto start2 = high_resolution_clock::now();
    
    vec eigval = eig_sym(A);
    
    auto finish2 = high_resolution_clock::now();
    
    std::cout << "Eigenvalues armadillo:" << std::endl; 
    print_and_sort(5,eigval);
    std::cout << "time(eig_sym(armadillo)):" 
        << duration_cast<milliseconds>(finish2 - start2).count() << "ms" 
        << std::endl;
    std::cout << "\n";

    ///////////////////////////////////////////////////////
    // new form of potensial V = omega_r^2*rho^2 + 1/rho //
    ///////////////////////////////////////////////////////
    double omega_r = 0.01; //constant defined
    vec V2(n, fill::zeros); 
    for(int i=0; i<n; i++) {
        V2(i) = omega_r*omega_r*rho(i)*rho(i) + 1./rho(i);
    }

    //diagonal vector
    vec d2(n); d2.fill(2./(h*h)); 
    for(int i=0; i<n; i++) {
        d2(i) += V2(i);
    }    

    //create tridiagonal matrix A for eigenvalues
    mat A2(n,n, fill::zeros);
    A2.diag() = d2; A2.diag(-1) = e; A2.diag(1) = e; 

    //create matrix R for eigenvectors
    mat R2(n,n, fill::zeros); R2.diag().fill(1.);

    //using armadillo to calculate eigenvalues
    vec eigval2; mat eigvecs2; 
    eig_sym(eigval2, eigvecs2, A2);
    vec eigvec21 = eigvecs2.col(0);
    vec eigvec22 = eigvecs2.col(1);
    vec eigvec23 = eigvecs2.col(2);
    std::cout << "Eigenvalues Oscillating:" << std::endl;
    print_and_sort(5,eigval2);
    
    //const char* filename = argv[2];
    //write_values_to_file(filename, eigval2, eigvec21, eigvec22, eigvec23, 
            //rho_max, rho_min, rho, n);
    
    return 0;
}
