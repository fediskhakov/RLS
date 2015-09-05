/*        John Rust, University of Maryland
 *        Bertel Schjerning, University of Copenhagen
 *        Fedor Iskhakov, University Technology Sidney
 *        February 2011
 */
#include "stdio.h"
#include <stdlib.h>
#include "math.h"
#include "mex.h"
#include "matrix.h"
//#include <omp.h>
#define MAX(A,B) ((A)>(B)?(A):(B))
#define MIN(A,B) ((A)>(B)?(B):(A))

// Set equal to the number of cores on your computer
#define NUM_THREADS 4

// Maximum number of equilibria for one combination of c1,c2
#define  MAXEQB 3
#define  MAXEQBINFO (MAXEQB+2)
#define nF 2

#define SYMMETRYONEDGESoff

#define MAXP 10
#define PI 3.14159265358979323846264338327950288419716939


/********************************************************/
// Declare data types
/********************************************************/

typedef struct bnestruct{
    double* ic1; double* ic2;
    double* c1; double* c2;
    double* p1; double* p2;
    double* pf1; double* pf2;
    double* s1; double* s2;
} Bnestruct; //structure of the pointers to output table columns

/*
typedef struct brstruct{
    double* ic1; double* ic2;
    double* c1; double* c2;
    double* A10; double* A20;
    double* A11; double* A21;
    double* B1; 
    double* C10; double* C20;
    double* C11; double* C21;
    double* a10; double* a20;
    double* a11; double* a21;
    double* a12; double* a22;
} Brstruct; //structure of the pointers to output table columns
*/

typedef struct brstruct{
    double* ic1; double* ic2; //indeces for the c1 and c2
    double* c1; double* c2;   //values of c1 and c2
    double* A0[2];      //coef for the value of not investing (alfa1) for firm 1 and 2
    double* A1[2];      //coef for the value of not investing (alfa2) for firm 1 and 2
    double* B1;         //df*(1-pfi) 
    double* C0[2];      //coef for the value of investing (gamma0) for firm 1 and 2          
    double* C1[2];      //coef for the value of investing (gamma1) for firm 1 and 2
    double* D0[2];      //coef to polinomial with P^0 for firm 1 and 2
    double* D1[2];      //coef to polinomial with P^1 for firm 1 and 2
    double* D2[2];      //coef to polinomial with P^2 for firm 1 and 2
} Brstruct; //coef for the analytical best responce and value functions

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

typedef struct mvstruct{
    double pf[2]; double kv; double logsumK; 
    double pti; double *h1; double *h2;
    double* eqbinfo; int neqbinfo;
} MVstruct; //structure of temp variables (for argument passing)



/********************************************************/
// Declare functions
// (so that they can appear in any order in the code)
/********************************************************/
void leap(MPstruct mp, Gamestruct* g, Bnestruct bne, Brstruct* br);//declaration
int id(int nC,int ic1,int ic2,int iC,int ieqb);
int ih(int nC, int ic1, int ic2, int iC);







//Proper index of a row in output matrix defined by ic1,ic2,ieqb,iC,nC
#define IdOUT ((ic1-iC)*(nC-iC)*MAXEQB+(ic2-iC)*MAXEQB+ieqb)
#define IdBNE (ic1*nC+ic2)
#define brOUT ((ic1-iC)*(nC-iC)+(ic2-iC))
#define NVarEqinfo 7
#define PRINTOUT ((ic1==2) && (ic2==2) && (iC==0))
int id(int nC,int ic1,int ic2,int iC,int ieqb) {
/*--------------------------------------------------------------------------------------
 Function id() return the row index in the output table given the indeces of costs
 and equilibrium.  See definition of output table "solution"
 **--------------------------------------------------------------------------------------*/
    return (ic1-iC)*(nC-iC)*MAXEQB+(ic2-iC)*MAXEQB+ieqb;
}
int ih(int nC, int ic1, int ic2, int iC) {
/*--------------------------------------------------------------------------------------
 Function ih() return the index of temp h array in games with c>0
**--------------------------------------------------------------------------------------*/
    return (ic1-iC)*(nC-iC)+(ic2-iC);
}

#include "functions.c"
#include "modelparts.c"
#include "analytical.c"





void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
     const char *keys[] = {"solution","c"};
     const int nkeys = 2;
     const char *brkeys[] = {"br"};
     const int nbrkeys = 1;
     Gamestruct *g;
     Bnestruct bne; 
     Brstruct *br; 
     MPstruct mp; 
     int j, iC, nC;
     
     if (nrhs != 2) mexErrMsgTxt("Error in seg.c: wrong number of inputs!");
     if (nlhs != 3) mexErrMsgTxt("Error in seg.c: wrong number of outputs!");
     if ((mxGetM(prhs[0])!=1) && (mxGetN(prhs[1])!=8)) mexErrMsgTxt("Wrong parameter vector, params!");
     if ((mxGetM(prhs[1])!=1) && (mxGetN(prhs[0])!=7)) mexErrMsgTxt("Wrong parameter vector, modparams!");
     
    /********************************************************/
    // UNZIP INPUTS
    /********************************************************/
    //   Parameters that index spaces, algoriths etc.
    //   (see setup for description of par)
    mp.nC       =   (int)   mxGetScalar(mxGetField(prhs[0],0,"nC")); 
    mp.nP       =   (int)   mxGetScalar(mxGetField(prhs[0],0,"nP")); 
    mp.og       =   (int)   mxGetScalar(mxGetField(prhs[0],0,"og"));
    mp.maxit    =   (int)   mxGetScalar(mxGetField(prhs[0],0,"maxit"));
    mp.ctol     = (double)  mxGetScalar(mxGetField(prhs[0],0,"ctol"));
    mp.cmin     = (double)  mxGetScalar(mxGetField(prhs[0],0,"cmin"));
    mp.cmax     = (double)  mxGetScalar(mxGetField(prhs[0],0,"cmax"));
    nC          =           mp.nC;   // very much used => for easy referencing
    //   model parameters
    //   (see setup.m for description mp)    
    mp.df       =   mxGetScalar(mxGetField(prhs[1],0,"df")); 
    mp.k1       =   mxGetScalar(mxGetField(prhs[1],0,"k1"));  
    mp.k2       =   mxGetScalar(mxGetField(prhs[1],0,"k2"));  
    mp.sigma    =   mxGetScalar(mxGetField(prhs[1],0,"sigma"));  
    mp.c_og     =   mxGetScalar(mxGetField(prhs[1],0,"c_og")); 
    mp.eta      =   mxGetScalar(mxGetField(prhs[1],0,"eta"));  
    mp.c_tr     =   mxGetScalar(mxGetField(prhs[1],0,"c_tr")); 
    //Grid for costs
    mp.cgrid= (double *) calloc(nC,    sizeof(double));     // grid of c on [cmin,cmax]
    for (j=0;j<nC;j++) mp.cgrid[j]=mp.cmin+(double)j*(mp.cmax-mp.cmin)/(nC-1);
         
    /********************************************************/
    // OUTPUT: bne structure
    /********************************************************/       
    // create nC^2 x 10  output Matlab matrix
    plhs[0]=mxCreateDoubleMatrix(nC*nC, 10, mxREAL);
    //Column names in output matrix
    //Use index ic1*nC+ic2 = IdBNE
    bne.ic1 =mxGetPr(plhs[0]);
    bne.ic2 =bne.ic1+nC*nC*1;
    bne.c1  =bne.ic1+nC*nC*2;
    bne.c2  =bne.ic1+nC*nC*3;
    bne.pf1 =bne.ic1+nC*nC*4;
    bne.pf2 =bne.ic1+nC*nC*5;
    bne.p1  =bne.ic1+nC*nC*6;
    bne.p2  =bne.ic1+nC*nC*7;
    bne.s1  =bne.ic1+nC*nC*8;
    bne.s2  =bne.ic1+nC*nC*9;
    //set output to NaN initially (use pointer to first element in output matrix)
    for (j=0;j<nC*nC*10;bne.ic1[j++]=mxGetNaN());
    //call main routine to fill up prepared output space
    s_bne(&mp, &bne);
    
    /********************************************************/
    // OUTPUT: br structure
    /********************************************************/       
    // structure 1xnC of nC^2 x 19  matrices
    plhs[1]=mxCreateStructMatrix(1, nC, nbrkeys, brkeys);
    br = calloc(nC, sizeof(Brstruct));
    //make fields
    for (iC=0;iC<nC;iC++) {
         mxSetField(plhs[1], iC, brkeys[0], mxCreateDoubleMatrix((nC-iC)*(nC-iC), 19, mxREAL));
        //Column names in output matrix
        //Use index (ic1-iC)*(nC-iC)+(ic2-iC) = brOUT
        //Below are pointer to different columts of the solution matrix to be output to Matlab
        br[iC].ic1=(double *) mxGetData(mxGetField(plhs[1], iC, brkeys[0]));
        br[iC].ic2=br[iC].ic1+(nC-iC)*(nC-iC)*1; 
        br[iC].c1=br[iC].ic1+(nC-iC)*(nC-iC)*2; 
        br[iC].c2=br[iC].ic1+(nC-iC)*(nC-iC)*3;
        br[iC].A0[0] =br[iC].ic1+(nC-iC)*(nC-iC)*4;
        br[iC].A0[1] =br[iC].ic1+(nC-iC)*(nC-iC)*5;
        br[iC].A1[0] =br[iC].ic1+(nC-iC)*(nC-iC)*6;
        br[iC].A1[1] =br[iC].ic1+(nC-iC)*(nC-iC)*7;
        br[iC].B1    =br[iC].ic1 +(nC-iC)*(nC-iC)*8;
        br[iC].C0[0] =br[iC].ic1+(nC-iC)*(nC-iC)*9;
        br[iC].C0[1] =br[iC].ic1+(nC-iC)*(nC-iC)*10;
        br[iC].C1[0] =br[iC].ic1+(nC-iC)*(nC-iC)*11;
        br[iC].C1[1] =br[iC].ic1+(nC-iC)*(nC-iC)*12;
        br[iC].D0[0] =br[iC].ic1+(nC-iC)*(nC-iC)*13;
        br[iC].D0[1] =br[iC].ic1+(nC-iC)*(nC-iC)*14;
        br[iC].D1[0] =br[iC].ic1+(nC-iC)*(nC-iC)*15;
        br[iC].D1[1] =br[iC].ic1+(nC-iC)*(nC-iC)*16;
        br[iC].D2[0] =br[iC].ic1+(nC-iC)*(nC-iC)*17;
        br[iC].D2[1] =br[iC].ic1+(nC-iC)*(nC-iC)*18;
        //set output to NaN initially (use pointer to first element in output matrix
        for (j=0;j<(nC-iC)*(nC-iC)*19;br[iC].ic1[j++]=mxGetNaN());
    }
    
    /********************************************************/
    // OUTPUT: g structure
    /********************************************************/
    plhs[2]=mxCreateStructMatrix(1, nC, nkeys, keys);
    g = calloc(nC, sizeof(Gamestruct));
     //create fields in output Matlab structure
     for (iC=0;iC<nC;iC++) {
         //i refers to different values of c (state of the art cost)
         mxSetField(plhs[2], iC, keys[0], mxCreateDoubleMatrix((nC-iC)*(nC-iC)*MAXEQB, 15, mxREAL));
         mxSetField(plhs[2], iC, keys[1], mxCreateDoubleScalar(mxGetNaN()));
         
         //Column names in output matrix solution
         //Use index ((ic1-iC)*(nC-iC)*MAXEQB+(ic2-iC)*MAXEQB+ieqb)
         //See function id(*,*,*,*) and #define IdOUT
         //Below are pointer to different columts of the solution matrix to be output to Matlab
         g[iC].ic1=(double *) mxGetData(mxGetField(plhs[2], iC, keys[0])); //key= solution
         g[iC].ic2=g[iC].ic1+(nC-iC)*(nC-iC)*MAXEQB;
         g[iC].ieqb=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*2;
         g[iC].c1=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*3;
         g[iC].c2=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*4;
         g[iC].eqbtype=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*5;
         /* Equilibrium type
          * 0= pure stategy
          * 1= mixed strategy */
         g[iC].v10=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*6;
         g[iC].v11=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*7;
         g[iC].v20=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*8;
         g[iC].v21=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*9;
         g[iC].p1=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*10;
         g[iC].p2=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*11;
         g[iC].seleqb=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*12;
         g[iC].pf1=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*13;
         g[iC].pf2=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*14;
         
         g[iC].c=(double *) mxGetData(mxGetField(plhs[2], iC, keys[1])); //key= c
         
         //set output to NaN initially (use pointer to first element in output matrix
         for (j=0;j<(nC-iC)*(nC-iC)*MAXEQB*15;g[iC].ic1[j++]=mxGetNaN());
     }

     //call main routine to fill up prepared output space
       leap(mp, g, bne, br);

     //cleanup
     free(mp.cgrid);
     free(g);
     free(br);
}












void leap(MPstruct mp, Gamestruct* g, Bnestruct bne, Brstruct* br) {
/********************************************************/
// SECTION 1: DECLARATIONS
/********************************************************/
    // Params and Modparams
    int nC; //Main parameter
    MVstruct mv;
    
    //double v, v1, r, p1;
    int i, j, iC, ic1, ic2, ieqb, neqb;
    //double q; //temps for c>0
    //double v, v1, r, p1;
    //int i, j, iF, iI, iC, ic1, ic2, ieqb, it, ip, neqb;
    double v, q;


    nC        =  mp.nC;               // very much used => for easy referencing
       
/********************************************************/
// SECTION 2: MEMORY ALLOCATION
/********************************************************/
    mv.h1      = (double *) calloc(nC*nC, sizeof(double));     // temp storage of formula parts, precomputing, see c>0
    mv.h2      = (double *) calloc(nC*nC, sizeof(double));     // temp storage of formula parts, precomputing, see c>0
    mv.eqbinfo = (double *) calloc(5*MAXEQBINFO, sizeof(double));     // storage for found equilibria in findeqb



/********************************************************/
// SECTION 3: LOOP OVER iC
/********************************************************/
for (iC=0;iC<nC;iC++) {
//for (iC=0;iC<1;iC++) {

    /********************************************************/
    // Initialize for iC 
    /********************************************************/
    g[iC].c[0]=mp.cgrid[iC]; // save value of c in output
    mv.pti=f_pti(iC,&mp);    // probability of a technological improvement. Initially improvements are assumed to be 
                             // 'incremental', that is the jump, if it occurs, is to c'=cgrid(i-1) with probability pti
                             // for iC=0 pti=0!
    mv.kv=f_kf(iC, &mp);     // Investment costs
    mv.logsumK=f_lnsum(2, mp.eta, 0, -mv.kv);        // Update mv.kv using investments costs, K(c)

    // precompute h(c1,c2,c) first element of H(c1,c2,c);
    if (iC==0) for (i=0;i<nC*nC;i++) {mv.h1[i]=0;mv.h2[i]=0;} 
    else {
//        omp_set_num_threads(NUM_THREADS);
//        #pragma omp parallel shared(mv, g, mp,iC,  nC) private(ic1, ic2, neqb, i, j) 
        {
//               printf("There are %d threads\n",omp_get_num_threads());
//        #pragma omp for
            for (ic1=iC;ic1<nC;ic1++) {
                for (ic2=iC;ic2<nC;ic2++) {
                    neqb=f_neqb(nC, ic1, ic2, iC-1, g);    
                    i=id(nC, ic1, ic2, iC-1, f_ers(iC-1, ic1, ic2, neqb));    // index for selected equilibrium
                    j=(ic1-iC)*(nC-iC)+(ic2-iC);    // index used for storage in h
                    mv.h1[j]=f_lnsum(nF, mp.eta, g[iC-1].v10[i], g[iC-1].v11[i]);
                    mv.h2[j]=f_lnsum(nF, mp.eta, g[iC-1].v20[i], g[iC-1].v21[i]);
                }
            }
        }
    }
    /********************************************************/
    // SECTION 3.0: Solve the (c,c,c) corner game
    /********************************************************/
    ic1=iC; ic2=iC; ieqb=0; // Start c=c1=c2=0, only one equilibrium in the corner  endgame

    // Solve (c,c,c) endgame for firm 1
    v=(bne.pf1[IdBNE]+mp.df*mv.pti*mv.h1[ih(nC, iC, iC, iC)]+mp.df*(1-mv.pti)*mv.logsumK)/(1-mp.df*(1-mv.pti));    
    // Save output
    g[iC].pf1[IdOUT]=bne.pf1[IdBNE];
    g[iC].pf2[IdOUT]=bne.pf2[IdBNE];
    g[iC].ic1[IdOUT]=ic1;
    g[iC].ic2[IdOUT]=ic2;
    g[iC].ieqb[IdOUT]=ieqb;
    g[iC].c1[IdOUT]=mp.cgrid[ic1];
    g[iC].c2[IdOUT]=mp.cgrid[ic2];
    g[iC].eqbtype[IdOUT]=1;    // 0 = pure strategy equilibrium
    g[iC].v10[IdOUT]=v;
    g[iC].v11[IdOUT]=v-mv.kv;
    g[iC].v20[IdOUT]=v;
    g[iC].v21[IdOUT]=v-mv.kv;
    g[iC].p1[IdOUT]=f_logit(nF, mp.eta, NULL, g[iC].v11[IdOUT], g[iC].v10[IdOUT]);
    g[iC].p2[IdOUT]=g[iC].p1[IdOUT];
    g[iC].seleqb[id(nC,ic1,ic2,iC,ieqb)]=1;    
    
    /**************************************************************************/
    // SECTION 3.1: solve the (c1,c,c) edge game
    /**************************************************************************/
    ic2=iC; ieqb=0; //only one equilibrium in the edge  endgame
    for (ic1=iC+1;ic1<nC;ic1++) {
        // Solve for BN-equaliprium prices and profits
        // Part of formula ??
        q=mv.pti*mv.h1[ih(nC, iC, iC, iC)]+(1-mv.pti)*g[iC].v10[id(nC, iC, iC, iC, 0)];
        // Save output
        g[iC].pf1[IdOUT]=bne.pf1[IdBNE];
        g[iC].pf2[IdOUT]=bne.pf2[IdBNE];
        g[iC].ic1[IdOUT]=ic1;
        g[iC].ic2[IdOUT]=ic2;
        g[iC].ieqb[IdOUT]=ieqb;
        g[iC].c1[IdOUT]=mp.cgrid[ic1];
        g[iC].c2[IdOUT]=mp.cgrid[ic2];
        g[iC].eqbtype[IdOUT]=1;    //mixed strategy equilibrium
        g[iC].v11[IdOUT]=bne.pf1[IdBNE]-mv.kv+mp.df*q+mp.df*(1-mv.pti)*mv.logsumK;
        g[iC].v10[IdOUT]=f_solvevf(
                bne.pf1[IdBNE]+mp.df*mv.pti*mv.h1[ih(nC, ic1, iC, iC)], 
                (1-mv.pti)*mp.df, 
                g[iC].v11[id(nC, ic1, iC, iC, 0)], 
                mp.eta, &mp);  
        g[iC].p1[IdOUT]=f_logit(nF, mp.eta, NULL, g[iC].v11[IdOUT], g[iC].v10[IdOUT]);
        g[iC].seleqb[IdOUT]=1;
#ifdef SYMMETRYONEDGES
        //by symmetry
        g[iC].v20[id(nC, ic2, ic1, iC, 0)]=g[iC].v10[IdOUT];
        g[iC].v21[id(nC, ic2, ic1, iC, 0)]=g[iC].v11[IdOUT];
        g[iC].p2[id(nC, ic2, ic1, iC, 0)]=g[iC].p1[IdOUT];
#else
        //no symmetry
        q=mv.pti*mv.h2[ih(nC, iC, iC, iC)]+(1-mv.pti)*g[iC].v20[id(nC, iC, iC, iC, 0)];
        g[iC].v20[IdOUT]=(bne.pf2[IdBNE]
                +mp.df*g[iC].p1[IdOUT]*q
                +mp.df*(1-g[iC].p1[IdOUT])*mv.pti*mv.h2[ih(nC, ic1, iC, iC)]
                +mp.df*(1-mv.pti)*mv.logsumK
                         )/(1-(1-mv.pti)*mp.df*(1-g[iC].p1[IdOUT]));   // equation (26)
        g[iC].v21[IdOUT]=g[iC].v20[IdOUT]-mv.kv;
        g[iC].p2[IdOUT]=f_logit(nF, mp.eta, NULL, -mv.kv, 0);;
#endif        
     }
    
    /**************************************************************************/
    // SECTION 3.2: solve the (c,c2,c) edge game
    /**************************************************************************/
    ic1=iC; ieqb=0; //only one equilibrium in the edge endgame
    for (ic2=iC+1;ic2<nC;ic2++) {
#ifndef SYMMETRYONEDGES
        //no symmetry
        q=mv.pti*mv.h2[ih(nC, iC, iC, iC)]+(1-mv.pti)*g[iC].v20[id(nC, iC, iC, iC, 0)];
        g[iC].v21[IdOUT]=bne.pf2[IdBNE]-mv.kv+mp.df*q+mp.df*(1-mv.pti)*mv.logsumK;
        g[iC].v20[IdOUT]=f_solvevf(
                bne.pf2[IdBNE]+mp.df*mv.pti*mv.h2[ih(nC, iC, ic2, iC)], 
                (1-mv.pti)*mp.df, 
                g[iC].v21[IdOUT], 
                mp.eta, &mp);
        g[iC].p2[IdOUT]=f_logit(nF, mp.eta, NULL, g[iC].v21[IdOUT], g[iC].v20[IdOUT]);
#endif
        q=mv.pti*mv.h1[ih(nC, iC, iC, iC)]+(1-mv.pti)*g[iC].v10[id(nC, iC, iC, iC, 0)];
        // Save output
        g[iC].pf1[IdOUT]=bne.pf1[IdBNE];
        g[iC].pf2[IdOUT]=bne.pf2[IdBNE];
        g[iC].ic1[IdOUT]=ic1;
        g[iC].ic2[IdOUT]=ic2;
        g[iC].ieqb[IdOUT]=ieqb;
        g[iC].c1[IdOUT]=mp.cgrid[ic1];
        g[iC].c2[IdOUT]=mp.cgrid[ic2];
        g[iC].eqbtype[IdOUT]=1;    //mixed strategy equilibrium
        //already have probabilities p2 of investment by firm 2
        g[iC].v10[IdOUT]=(
                bne.pf1[IdBNE]
                +mp.df*g[iC].p2[id(nC,iC,ic2,iC,0)]*q
                +mp.df*(1-g[iC].p2[id(nC,iC,ic2,iC,0)])*mv.pti*mv.h1[ih(nC, iC, ic2, iC)]
                +mp.df*(1-mv.pti)*mv.logsumK
                )/(1-(1-mv.pti)*mp.df*(1-g[iC].p2[id(nC,iC,ic2,iC,0)]));   // equation (26)
        g[iC].v11[IdOUT]=g[iC].v10[IdOUT]-mv.kv;
        g[iC].p1[IdOUT]=f_logit(nF, mp.eta, NULL, g[iC].v11[IdOUT], g[iC].v10[IdOUT]);
        g[iC].seleqb[IdOUT]=1;
#ifdef SYMMETRYONEDGES
        //by symmetry
        g[iC].v20[id(nC, ic2, ic1, iC, 0)]=g[iC].v10[IdOUT];
        g[iC].v21[id(nC, ic2, ic1, iC, 0)]=g[iC].v11[IdOUT];
        g[iC].p2[id(nC, ic2, ic1, iC, 0)]=g[iC].p1[IdOUT];
#endif
    }
                  


    /**************************************************************************/
    // SECTION 3.3: solve the (c1,c2,c) interior game
    /**************************************************************************/
//     omp_set_num_threads(NUM_THREADS);
//     #pragma omp parallel shared(mv, g, mp, nC, iC, bne) private(ic1, ic2, i1, i2, A0, A1, B1, C0, C1, a0, a1, a2, pstar, nroot, ieqb, i, j, A, B, v0_0, v0_1, pbr, save, neqb, nmeqb) 
//     {
// //     printf("There are %d threads\n",omp_get_num_threads());
//     #pragma omp for 

    //loop over the interier points
    for (ic1=iC+1;ic1<nC;ic1++) {
        for (ic2=iC+1;ic2<nC;ic2++) {
            //calculate analytical solutions and save to proper output structures
            analytical_solutions (ic1,ic2,iC,mp,mv,bne,g,br);
        }
    }
//     }  // end parallel execution
} //continue loop over iC
    
    free(mv.h1);
    free(mv.h2);
    free(mv.eqbinfo);
}
