#define MAX(A,B) ((A)>(B)?(A):(B))
#define MIN(A,B) ((A)>(B)?(B):(A))

#include "mex.h"
#include "math.h"
#include "stdio.h"


typedef struct mpstruct{
    double eta; double sigma;
    double df; double k1; double k2;  double c_og; double c_tr;
    int nC; double cmin; double cmax; double *cgrid;
    int og;
    int maxit; double tolerance; int nP;
} MPstruct; //structure of model parameters (for argument passing)
typedef struct mvstruct{
    double pf[2]; double kv; double logsumK; double logsum1; double logsum2; 
    double pti; double *h1; double *h2;
} MVstruct; //structure of temp variables (for argument passing)

#include "functions.c"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double eta;
    double *x;
    double *y;
    int n;
    
    
     if (nrhs != 2) mexErrMsgTxt("Error in lnsum.c: wrong number of inputs!");
     if (nlhs != 1) mexErrMsgTxt("Error in lnsum.c: wrong number of outputs!");
     if (mxGetN(prhs[1])!=1) mexErrMsgTxt("Error in lnsum.c: x must be column vector!");
    
     //input:
     n= (int) mxGetM(prhs[0]); 
     if (n > 10) mexErrMsgTxt("lnsum cannot be evaluated for n>10");

     x=mxGetPr(prhs[0]);
     eta=mxGetScalar(prhs[1]);

     //output
     plhs[0]=mxCreateDoubleScalar(mxREAL);
     y=mxGetPr(plhs[0]);

     if      (n==1) y[0]=f_lnsum(n, eta, x[0]);
     else if (n==2) y[0]=f_lnsum(n, eta, x[0], x[1]);
     else if (n==3) y[0]=f_lnsum(n, eta, x[0], x[1], x[2]);
     else if (n==4) y[0]=f_lnsum(n, eta, x[0], x[1], x[2], x[3]);
     else if (n==5) y[0]=f_lnsum(n, eta, x[0], x[1], x[2], x[3], x[4]);
     else if (n==6) y[0]=f_lnsum(n, eta, x[0], x[1], x[2], x[3], x[4], x[5]);
     else if (n==7) y[0]=f_lnsum(n, eta, x[0], x[1], x[2], x[3], x[4], x[5], x[6]);
     else if (n==8) y[0]=f_lnsum(n, eta, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]);
     else if (n==9) y[0]=f_lnsum(n, eta, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]);
     else if (n==10) y[0]=f_lnsum(n, eta, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]);
    
}
