/* Outputting to seperate file and plotting in python */
#include <cmath> //exp, log10
#include <armadillo> //vec, mat, solve, lu
#include <string> //atoi
#include <cstdlib> //system
#include <iostream> //cout
#include <chrono> //high_resolution_clock

using namespace arma;
using namespace std::chrono;

vec solver(int n, vec a, vec b, vec c, vec bfu) {
  /* solving with tridiagonal algorithm */

  //initialise vectors
  vec v(n+2, fill::zeros); 
  vec temp(n+2,fill::zeros);
  vec btemp(n+2, fill::zeros);
  
  //start timer(from std::chrono)
  auto start = high_resolution_clock::now();
  
  //the algorithm
  btemp(1) = b(1);
  v(1) = bfu(1)/btemp(1);
  temp(1) = c(1)/btemp(1);
  for(int i=2; i<n+1 ;i++) {
    btemp(i) = b(i) - a(i)*temp(i-1);
    temp(i) = c(i)/btemp(i);
    v(i) = (bfu(i) - a(i)*v(i-1))/btemp(i);
  }
  
  for(int i=n; i>0; i--) {
    v(i) -= temp(i)*v(i+1);
  }
  
  //stop timer(from std::chrono)
  auto finish = high_resolution_clock::now();
  
  //print timer results
  std::cout << "Tridiagonal:" 
	    << duration_cast<nanoseconds>(finish - start).count()
	    << "ns" << '\n';
  
  return v;
}

vec LU_solve(int n) {
  /* solve the problem with LU-decomposition */
  double h = 1. / (n-1); //step-length
  
  //matrix as given in project 
  mat A(n,n, fill::zeros); 
  for(int i=0; i<n-1 ;i++) {
    /* set values in matrix */
    A(i,i+1) = -1; A(i,i) = 2; A(i+1,i) = -1;
  }
  
  A(n-1,n-1) = 2; //set last value(in vector b)
  
  //right-side of system(b-tiled)
  vec bfu(n, fill::zeros);
  for(int i=0; i<n-1; i++) bfu[i] = h*h*100*exp(-10*i*h);
  
  mat L, U, P; //make matrix(upper, lower and reduced)

  //start clock(LU)
  auto start2 = high_resolution_clock::now();
 
  lu(L,U,P,A); //lu-decomposition in armadillo
  
  //solve LU with armadillo
  vec LU_y = P.t()*solve(L,bfu); vec LU_x = solve(U,LU_y);
  
  //finish clock(LU)
  auto finish2 = high_resolution_clock::now();
  cout << "LU solve:" 
       << duration_cast<nanoseconds>(finish2 - start2).count()
       << "ns" << '\n';

  return LU_x;
}

vec bfunc(int n) {
  /* function for right side of system*/
  
  //set vectors
  vec bfuu(n+2, fill::zeros); 
  vec xbf(n+2, fill::zeros);
 
  double h = 1./(n+1); //step-length
  
  //initial values
  xbf(1) = 1*h; bfuu(1) = h*h*100*exp(-10*xbf(1)); 
  for(int i=2; i<n+1 ;i++) {
    /* assign values to x and bfunc */
    xbf(i) = i*h;
    bfuu(i) = h*h*100*std::exp(-10*xbf(i));
  }
  xbf(n+1) = 1;
  return bfuu;
}
double relerr(int n, vec vres) {
  /* function for solving relative error */
  vec u(n+2, fill::zeros); //vector for exact solution 
  double h = 1./(n+1); //step-lenght
  
  u(1) = 1 - (1 - exp(-10))*1*h - exp(-10*1*h);
  for(int i=2; i<n+1 ;i++) {
    /* assign values to u */
    u(i) = 1 - (1 - exp(-10))*i*h - exp(-10*i*h);
  }

  double epsilon = -1; double eps_temp = -1;
  for(int i=1; i<n-1 ;i++) {
    /* assign values to epsilon */
    eps_temp = std::abs((vres(i) - u(i))/u(i));
    if(eps_temp > epsilon) epsilon = eps_temp;
  }
  
  return std::log10(epsilon);
}

void opener(const char* filename, vec res, vec x, int n, bool err) {
  /* function for writing to file */

  std::ofstream outputfile; //create filetype
  outputfile.open(filename); //create(and/or open) file
  
  if(err == false) {
    //x-values for exact solution(for right length of vector)
    int t = n+2;
    outputfile << "n = " << t << '\n' << res << x; //output to file
    outputfile.close(); //close file
  }
  else {
    outputfile << res;
    outputfile.close();
  }
}

int main(int argc, char* argv[]) {
  if(argc < 3) {
    cout << "USAGE: mode n, mode: 1=tri, 2=lu, 3=error" << std::endl;
  }
  
  int n = std::atoi(argv[2]); //number of iterations(from cmdl)
  int mode = std::atoi(argv[1]); //mode(tri,lu,err)

  //set vectors in matrix
  vec a(n+2); a.fill(-1); a(0) = 0;
  vec b(n+2); b.fill(2);
  vec c(n+2); c.fill(-1); c(n+1) = 0;

  double h = 1./(n+1); //step-length

  vec x(n+2, fill::zeros); x(1) = 1*h; //x-value 
  for(int i=2; i<n+1 ;i++) x(i) = i*h; //set x-values
  x(n+1) = 1; //end boundary
  
  if(mode==1) { 
    /* run tridiagonal */
    vec res = solver(n, a, b, c, bfunc(n)); //get result
    opener("tri_res.dat", res, x, n, false); //make file
    std::system("python 2fe_plot.py tri_res.dat"); //run python plot
  }
  
  else if(mode==2) {
    /* run lu */
    vec lu_res = LU_solve(n);
  }
  
  else if(mode==3) {
    /* run error */
    const int dim = 6; //dimension for n-values
    int N[dim] = {10,100,1000,10000,100000,1000000};
    
    vec eps(dim, fill::zeros); //set eps
    for(int i=0; i<dim ;i++) {
      /* run through all n-values */
      /* this is really bad -.- */
      vec a_temp(N[i]+2); a_temp.fill(-1); a_temp[0] = 0;
      vec b_temp(N[i]+2); b_temp.fill(2);
      vec c_temp(N[i]+2); c_temp.fill(-1); c_temp[N[i]+1] = 0;
      vec bfu_temp = bfunc(N[i]);
      vec solve_temp = solver(N[i], a_temp, b_temp, c_temp, bfu_temp);
      eps[i] = relerr(N[i], solve_temp);
    }
    opener("tri_error.dat", eps, x, 1, true);
    std::system("python 2fe_plot.py tri_error.dat");
  }
  
  else {
    /* Standard output */
    cout << "specify" << std::endl;
  }

  return 0;
}
