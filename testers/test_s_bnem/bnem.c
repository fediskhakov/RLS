/*        John Rust, University of Maryland
 *        Bertel Schjerning, University of Copenhagen
 *        Fedor Iskhakov, University Technology Sidney
 *        February 2011
 */
#include "stdio.h"
#include "math.h"
#include "mex.h"
#include "matrix.h"
#include <omp.h>

#define MAX(A, B) ((A)>(B)?(A):(B))
#define MIN(A, B) ((A)>(B)?(B):(A))

// Set equal to the number of cores on your computer
#define NUM_THREADS 4  

// Maximum number of equilibria for one combination of c1,c2
#define  MAXEQB 3
#define nF 2
//Proper index of a row in output matrix defined by ic1,ic2,ieqb,iC,nC
#define IdOUT ((ic1-iC)*(nC-iC)*MAXEQB+(ic2-iC)*MAXEQB+ieqb)

typedef struct bnestruct{
    double* ic1; double* ic2;
    double* c1; double* c2;
    double* p1; double* p2;
    double* pf1; double* pf2;
    double* s1; double* s2;
} Bnestruct; //structure of the pointers to output table columns

typedef struct mpstruct{
    double eta; double sigma;
    double df; double k1; double k2;  double c_og; double c_tr;
    int nC; double cmin;  double cmax; double *cgrid;
    int og;
    int maxit; double ctol; int nP;
} MPstruct; //structure of model parameters (for argument passing)

typedef struct mvstruct{
    double pf[2];//space for firm profits
    double s[3]; //space for market shares for nF+og=2+1
    double p[3]; //vector of prices
    double kv; double logsumK; double logsum1; double logsum2;
    double pti; double *h1; double *h2;
} MVstruct; //structure of temp variables (for argument passing)


int id(int nC, int ic1, int ic2, int iC, int ieqb) {
    /*--------------------------------------------------------------------------------------
     * Function id() return the row index in the output table given the indeces of costs
     * and equilibrium.  See definition of output table "solution"
     **--------------------------------------------------------------------------------------*/
    return (ic1-iC)*(nC-iC)*MAXEQB+(ic2-iC)*MAXEQB+ieqb;
}
int ih(int nC, int ic1, int ic2, int iC) {
    /*--------------------------------------------------------------------------------------
     * Function ih() return the index of temp h array in games with c>0
     **--------------------------------------------------------------------------------------*/
    return (ic1-iC)*(nC-iC)+(ic2-iC);
}

#include "functions.c"


void s_bne(MPstruct *mp, Bnestruct bne);//declaration

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    Bnestruct bne;
    int j, nC;
    
    MPstruct mp; // Param structures
    
    if (nrhs != 2) mexErrMsgTxt("Error in seg.c: wrong number of inputs!");
    if (nlhs != 1) mexErrMsgTxt("Error in seg.c: wrong number of outputs!");
    if ((mxGetM(prhs[0])!=1) && (mxGetN(prhs[1])!=8)) mexErrMsgTxt("Wrong parameter vector, params!");
    if ((mxGetM(prhs[1])!=1) && (mxGetN(prhs[0])!=7)) mexErrMsgTxt("Wrong parameter vector, modparams!");
    
    
    /********************************************************/
// UNZIP INPUTS
    /********************************************************/
//   Parameters that index spaces, algoriths etc.
//   (see setup for description of par)
    
    mp.nC       =   (int)   mxGetScalar(mxGetField(prhs[0], 0, "nC"));
    mp.nP       =   (int)   mxGetScalar(mxGetField(prhs[0], 0, "nP"));
    mp.og       =   (int)   mxGetScalar(mxGetField(prhs[0], 0, "og"));
    mp.maxit    =   (int)   mxGetScalar(mxGetField(prhs[0], 0, "maxit"));
    mp.ctol     =           mxGetScalar(mxGetField(prhs[0], 0, "ctol"));
    mp.cmin     =           mxGetScalar(mxGetField(prhs[0], 0, "cmin"));
    mp.cmax     =           mxGetScalar(mxGetField(prhs[0], 0, "cmax"));
    nC          =           mp.nC;   // very much used => for easy referencing
    
//   model parameters
//   (see setup for description mp)
    mp.df       =   mxGetScalar(mxGetField(prhs[1], 0, "df"));
    mp.k1       =   mxGetScalar(mxGetField(prhs[1], 0, "k1"));
    mp.k2       =   mxGetScalar(mxGetField(prhs[1], 0, "k2"));
    mp.sigma    =   mxGetScalar(mxGetField(prhs[1], 0, "sigma"));
    mp.c_og     =   mxGetScalar(mxGetField(prhs[1], 0, "c_og"));
    mp.eta      =   mxGetScalar(mxGetField(prhs[1], 0, "eta"));
    mp.c_tr     =   mxGetScalar(mxGetField(prhs[1], 0, "c_tr"));
    
    /********************************************************/
// OUTPUT: bne matrix
    /********************************************************/
    
    //C representation of Matlab output (pointers)
    
    // create nC^2 x 10  output Matlab matrix
    plhs[0]=mxCreateDoubleMatrix(nC*nC, 10, mxREAL);
    
    //Column names in output matrix
    //Use index ic1*nC+ic2
    //Below are pointer to different columts of the solution matrix to be output to Matlab
    bne.ic1=mxGetPr(plhs[0]);
    bne.ic2=bne.ic1+nC*nC*1;
    bne.c1=bne.ic1+nC*nC*2;
    bne.c2=bne.ic1+nC*nC*3;
    bne.pf1=bne.ic1+nC*nC*4;
    bne.pf2=bne.ic1+nC*nC*5;
    bne.p1 =bne.ic1+nC*nC*6;
    bne.p2 =bne.ic1+nC*nC*7;
    bne.s1 =bne.ic1+nC*nC*8;
    bne.s2 =bne.ic1+nC*nC*9;
    
    //set output to NaN initially (use pointer to first element in output matrix
    for (j=0;j<nC*nC*8;bne.ic1[j++]=mxGetNaN());
    
    //call main routine to fill up prepared output space
    s_bne(&mp, bne);
//       free(bne);
    
}

void s_bne(MPstruct *mp, Bnestruct bne) {
    /*--------------------------------------------------------------------------------------
     *        solves the "end game" equilibria when costs using the state of the
     *        art production technology have reached the lowest possible level, which is normalized
     *        to zero and are treated as an absorbing state
     **--------------------------------------------------------------------------------------*/
    
    /********************************************************/
    // SECTION 1: DECLARATIONS
    /********************************************************/
    // Params and Modparams
    int i, j, ic1, ic2;
    double cp, p[3], s[2], pf[2], c1, c2;
    int nC, it;
    double det, a[4], b[2];
    
    nC        =  mp[0].nC;               // very much used => for easy referencing
    
    
    /********************************************************/
    // SECTION 2: MEMORY ALLOCATION
    /********************************************************/
    mp[0].cgrid= (double *) calloc(nC,    sizeof(double));     // grid of c on [cmin,cmax]
    //Grid for costs
    for (i=0;i<nC;i++) mp[0].cgrid[i]=mp[0].cmin+(double)i*(mp[0].cmax-mp[0].cmin)/(nC-1);
    
    
    /********************************************************/
    // BN-equaliprium prices, profits and market shares
    /********************************************************/
    
        omp_set_num_threads(NUM_THREADS);
        #pragma omp parallel private(ic1, ic2, c1, c2, it, p,s,a,b,det, cp) shared(mp, bne) // shared(x, cols) reduction(+: sum) 
        {
  
    for (ic1=0;ic1<nC;ic1++) { //ic1
        #pragma omp for
//                      
        for (ic2=0;ic2<=ic1;ic2++) {  //ic2
            c1=mp[0].cgrid[ic1];    // very much used => for easy referencing
            c2=mp[0].cgrid[ic2];
            
            // Pure Bertrand solution
            f_logit(nF+mp[0].og, mp[0].sigma, s, (-1)*c1, (-1)*c2, (-1)*mp[0].c_og); //put logit probs into pi
            if (mp[0].og==0) { //no outside good
                p[0]=MAX(c1, c2);
            }
            else { //outside good
                if (c1>=c2 && c1>=mp[0].c_og) {
                    p[0]=MAX(c2, mp[0].c_og);
                }
                else if (c2>c1 && c2>mp[0].c_og) {
                    p[0]=MAX(c1, mp[0].c_og);
                }
                else if (mp[0].c_og>c1 && mp[0].c_og>c2) {
                    p[0]=MAX(c1, c2);
                }
            }
            p[1]=p[0];
            
            if (mp[0].sigma>0) { // if mp[0].sigma>0
                cp=2*mp[0].ctol;
                p[2]=p[0];
                /*find prices that satisfy F.O.C. of prifit maximization
                 * for both firms => solve 2x2 system of unlinear equations
                 * by Newton method:
                 * i=1..3: sigma-(1-s[i])(p[i]-c[i])=0, where s[i] is logit prob of full x
                 * Bertrand prices are used as startting values */
                
                for (it=0;((it<mp[0].maxit) && (cp>=mp[0].ctol));it++) { //iterations
                    f_logit(nF+mp[0].og, mp[0].sigma, s, (-1)*p[0], (-1)*p[1], (-1)*p[2]); //put logit probs into pi
                    //no outside good => (cross-)derivatives 2x2
                    a[0]=-(1-s[0])*(1+s[0]*(p[0]-c1)/mp[0].sigma);
                    a[1]=s[0]*s[1]*(p[0]-c1)/mp[0].sigma;
                    a[2]=s[0]*s[1]*(p[1]-c2)/mp[0].sigma;
                    a[3]=-(1-s[1])*(1+s[1]*(p[1]-c2)/mp[0].sigma);
                    det=a[0]*a[3]-a[1]*a[2];
                    //F.O.C.
                    b[0]=mp[0].sigma-(1-s[0])*(p[0]-c1);
                    b[1]=mp[0].sigma-(1-s[1])*(p[1]-c2);
                    p[0]-=  b[0]*a[3]/det - b[1]*a[1]/det;
                    p[1]-=- b[0]*a[2]/det + b[1]*a[0]/det;
                    
                    //criterion = MAX(ABS(F.O.C.))
                    cp=MAX(fabs(b[0]), fabs(b[1]));
                }
                if (it>=mp[0].maxit) mexErrMsgTxt("No convergence in s_bnem()!");
            }
            else if (mp[0].sigma==0) {
            }
            else mexErrMsgTxt("Don't know how to handle sigma<0!");
            pf[0]=(p[0]-c1)*s[0];
            pf[1]=(p[1]-c2)*s[1];
            
            // Save output
            i=ic1*nC+ic2;
            j=ic2*nC+ic1;
            
            bne.ic1[i]=ic1;
            bne.ic2[i]=ic2;
            bne.c1[i]=mp[0].cgrid[ic1];
            bne.c2[i]=mp[0].cgrid[ic2];
            
            bne.p1[i]=p[0];
            bne.p2[i]=p[1];
            bne.pf1[i]=pf[0];
            bne.pf2[i]=pf[1];
            bne.s1[i]=s[0];
            bne.s2[i]=s[1];
            
            // by symetry
            bne.ic1[j]=ic2;
            bne.ic2[j]=ic1;
            bne.c1[j]=mp[0].cgrid[ic2];
            bne.c2[j]=mp[0].cgrid[ic1];
            
            bne.p2[j]=bne.p1[i];
            bne.p1[j]=bne.p2[i];
            bne.pf2[j]=bne.pf1[i];
            bne.pf1[j]=bne.pf2[i];
            bne.s2[j]=bne.s1[i];
            bne.s1[j]=bne.s2[i];
        }
        }
    }
    free(mp[0].cgrid);
}