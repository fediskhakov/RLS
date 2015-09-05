/*        John Rust, University of Maryland
 *        Bertel Schjerning, University of Copenhagen
 *        Fedor Iskhakov, University Technology Sidney
 *        February 2011
 */
#include "stdio.h"
#include <stdlib.h>
#include "math.h"
#include "mex.h"

#define MAX(A, B) ((A)>(B)?(A):(B))
#define MIN(A, B) ((A)>(B)?(B):(A))

#define  MAXEQB 3
#define EQBMAX 10
#define MAXP 10
#define PI 3.14159265358979323846264338327950288419716939

typedef struct gamestruct{
    double* v10; double* v11; double* v20; double* v21; double* p1; double* p2;
    double* ic1; double* c1; double* ic2; double* c2; double* c;
    double* eqbtype; double* ieqb; double* seleqb; double* pf1; double* pf2;
} Gamestruct; //structure of the pointers to output table columns

typedef struct mpstruct{
    double eta; double sigma;
    double df; double k1; double k2;  double c_og; double c_tr;
    int nC; double cmin; double cmax; double *cgrid;
    int og;
    int maxit; double ctol; int nP;
} MPstruct; //structure of model parameters (for argument passing)

#include "analytical.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double* eqbinfo;
    if (nrhs != 0) mexErrMsgTxt("Error in seg.c: wrong number of inputs!");
    if (nlhs != 1) mexErrMsgTxt("Error in seg.c: wrong number of outputs!");
    
    /********************************************************/
    // OUTPUT: eqbinfo matrix
    /********************************************************/
    plhs[0]=mxCreateDoubleMatrix(5, EQBMAX, mxREAL);
    eqbinfo=mxGetPr(plhs[0]);
    
    findeqb(eqbinfo);
}

