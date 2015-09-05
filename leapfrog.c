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

int COUNT=0;

/********************************************************/
// IMPORTANT SWITCHES
/********************************************************/

// Set equal to the number of cores on your computer
#define NUM_THREADS 4  

// Maximum number of equilibria for one combination of c1,c2
#ifndef MAXEQB
    #define  MAXEQB 3
#endif
#define  MAXEQBINFO (MAXEQB+2)
#define nF 2

#define MAXP 10
#define PI 3.14159265358979323846264338327950288419716939
//tolerance used in output functions
#define OTOLERANCE 1e-8

#ifndef PRINTeqbloop
#define PRINTeqbloop 1
//PRINTeqbloop: 0 nothing
//              1 final info
//              2 iteration info
#endif
#ifndef PRINTeqbstr
#define PRINTeqbstr 0
//PRINTeqbstr: 0 no
//             1 all feasible eqbstrings
#endif

//calculate output statistics in eqbstraings cycle - takes more time
#define OUTSTATISTICS



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
    double* x10; double* x11; double* x20; double* x21; //additional value functions for alternating moves
    double* ic1; double* c1; double* ic2; double* c2; double* c;
    double* eqbtype; double* ieqb; double* seleqb; double* pf1; double* pf2;
    double* neqb;
    double* ecv10; double* ecv11; double* ecv20; double* ecv21; double* ecx10; double* ecx11; double* ecx20; double* ecx21; //additional values for expected cost
    double* ecm1, *ecm2; //total expected costs for both duopolists (in each point of state pyramid)
} Gamestruct; //structure of the pointers to output table columns

typedef struct mpstruct{
    double dt; double eta; double sigma;
    double df; double k1; double k2;  double c_og; double c_tr;
    double tpm11; double tpm12; double tpm21; double tpm22; // transistion probability for alternating moves
    int nC; double cmin; double cmax; double *cgrid;
    int og;
    int maxit; double ctol; int nP; double *pti; double *tpm;
    int esr; size_t esrmax; size_t esrmaxout; size_t esrstart; double * esrstartpt;
    int alternate; int analytical;
    double*out4;    //forth output - used in esr
    double *monsoluiton, *monvsown; //map of state space where monopoly invests
} MPstruct; //structure of model parameters (for argument passing)

typedef struct mvstruct{
    double pf[2]; double kv; double logsumK; double logsum1; double logsum2; 
    double pc; double *h1; double *h2;
    double* eqbinfo; int neqbinfo;
} MVstruct; //structure of temp variables (for argument passing)

typedef struct brc_amstruct{
    double pf[2]; double p[2]; double v0[2]; double v1[2]; double x0[2]; double x1[2]; 
    double f[2]; double A[2]; double B[2]; 
    double A0[2]; double B0[2]; 
    double A1[2]; double B1[2]; 
    double h11; double h12;  double h21; double h22;
    double h1; double h2;  
    double logsum1; double logsum2; double logsumK; double pc; double kv; 
    
    double betnpc;
} Brc_amstruct; //structure of temp variables for (for argument passing)


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
    return IdOUT;
          //(ic1-iC)*(nC-iC)*MAXEQB+(ic2-iC)*MAXEQB+ieqb;
}
int ih(int nC, int ic1, int ic2, int iC) {
/*--------------------------------------------------------------------------------------
 Function ih() return the index of temp h array in games with c>0
**--------------------------------------------------------------------------------------*/
    return brOUT;
    //(ic1-iC)*(nC-iC)+(ic2-iC);
}

//some function declarations (needed given the order of inclusion)
size_t lexindex1(int *eqstring,int nstr);
static void printeqb(int *eqstring,int neqstr,size_t indx,int *mask);

#include "functions.c"
#include "esr.c"
#include "modelparts.c"
#include "ffxp.c"
#include "brp_sm.c"
#include "analytical_sm.c"
#include "analytical_am.c"
#include "leap_sm.c"
#include "leap_am.c"


static void printeqb(int *eqstring,int neqstr,size_t indx,int *mask) {
    int j;
    int i,iC,ic1,ic2,nC=5;
    static first=0;
    if (first==0)
    {   //print ic,ic1,ic2
        printf("state space print out (infeasible points may be masked)\n");
        for (i=0;i<neqstr;i++) {iKinv (neqstr-1-i,nC,&iC,&ic1,&ic2);printf("%1d ",iC);}
        printf(" <-- iC\n");
        for (i=0;i<neqstr;i++) {iKinv (neqstr-1-i,nC,&iC,&ic1,&ic2);printf("%1d ",ic1);}
        printf(" <-- ic1\n");
        for (i=0;i<neqstr;i++) {iKinv (neqstr-1-i,nC,&iC,&ic1,&ic2);printf("%1d ",ic2);}
        printf(" <-- ic2\n");
        for (i=0;i<neqstr;i++) {printf("--");}
        printf("-- lexindex\n");
        first=1;
    }
    if (mask==NULL) for (j=0;j<neqstr;printf("%d ",eqstring[neqstr-1-j++]));
    else
    {
        for (j=0;j<neqstr;j++)
        {
            if (mask[neqstr-1-j]) printf("%d ",eqstring[neqstr-1-j]);
            else printf(". ");
        }
    }
    //printf conditional on compiler
#if defined(_MSC_VER)
    printf(" : %u\n",indx);
#elif defined(__GNUC__)
    printf(" : %zu\n",indx);
#endif
}

static int compout(const void *a,const void *b); //declare comparison function
int outputiter(size_t lexindex,int* eqstring,MPstruct *mp,Gamestruct *g,double *output) {
/*--------------------------------------------------------------------------------------
 Calculates all output variables
 Outputs the results from iteration to the output table
 Output is stored COLUMN BY COLUMN to allow for sorting
**--------------------------------------------------------------------------------------*/
    static size_t outindex=0; //the index of the column to output to next, KEEPS VALUE BETWEEN CALLS
    size_t outindex0;
    int length0;
    int i,nC,neqstr; //column number
    double *bsres;
    int iC,ic1,ic2,ieqb,jC,jp,j;
    int pure=0,symmetry=0,leapfrog=0,underinvestment=0;
    int *statepr, *stateprm; //marks of potential equilibrium path (pr>0 marked by 1)
    double p1,p2;

    //initialize call
    if (!eqstring)
    {   //initialize the output in the beginning of eqbstr cycle
        outindex=0;
        return 0;
    }   
    //initialize dimentions
    nC=mp[0].nC;
    neqstr=nC*(nC+1)*(2*nC+1)/6; //number of state points
    //first check if there is still space in output matrix
    if (outindex>mp[0].esrmaxout) return 2; //return err 2 qiuetly above the limit
    if (outindex==mp[0].esrmaxout) 
    {
        mexWarnMsgTxt("The number of distinct equilibria found is greater than sw.esrmax!");
        outindex++;
        return 1; //return err 1
    }
    
#ifdef OUTSTATISTICS
    //PART 1 Compute interesting statistics from g structure

    //Positive prob of state points
    statepr=calloc(neqstr,sizeof(int));
    stateprm=calloc(neqstr,sizeof(int));
    for (i=0;i<neqstr;statepr[i++]=0);
    for (i=0;i<neqstr;stateprm[i++]=0); //0=both firms may have the move, 1=only firm 1 may move here, 2=only firm 2 may move here
    statepr[iK(nC-1,nC-1,nC-1,nC)]=1;//assume that c0,c0,c0 is visited with certainty
    if (mp[0].alternate) stateprm[iK(nC-1,nC-1,nC-1,nC)]=1;//assume that at c0,c0,c0 is firm1 has the right to move
    //cycle over the state space from the top
    //printf("statepr pr>0 marks:\n");
    for (iC=nC-1;iC>=0;iC--)
    {
        // printeqb(statepr,neqstr,iC,NULL);
        // if (mp[0].alternate) printeqb(stateprm,neqstr,90+iC,NULL);
        if (iC<nC-1) //not very top layer
        {   //second to the top layer and lower
            //FIRST look at interior
            for (ic1=iC+1;ic1<nC;ic1++)
            {
                for (ic2=iC+1;ic2<nC;ic2++)
                {
                    jp=iK(iC,ic1,ic2,nC); //current point index
                    p1=g[iC].p1[id(nC,ic1,ic2,iC,eqstring[jp])];
                    p2=g[iC].p2[id(nC,ic1,ic2,iC,eqstring[jp])];
                    if (mp[0].alternate)
                    {   //alternating move
                        statepr[iK(iC,iC,ic2,nC)]=statepr[iK(iC,iC,ic2,nC)] || (statepr[jp] && p1>0 && stateprm[jp]!=2);
                        statepr[iK(iC,ic1,iC,nC)]=statepr[iK(iC,ic1,iC,nC)] || (statepr[jp] && p2>0 && stateprm[jp]!=1);
                        if (statepr[jp] && fabs(p1-1)<OTOLERANCE && stateprm[jp]!=2 && mp->pti[iC*nC+iC]<OTOLERANCE)
                        {   //only if technology certainly improves leave a trace of who had the move on this level
                            stateprm[iK(iC,iC,ic2,nC)]=stateprm[jp];
                        }
                        if (statepr[jp] && fabs(p1-1)<OTOLERANCE && stateprm[jp]!=2 && iC==0)
                        {   //only if technology certainly improves leave a trace of who had the move on this level
                            if (stateprm[jp]==1 && fabs(p1-1)<OTOLERANCE && fabs(mp->tpm11-1)<OTOLERANCE) stateprm[iK(iC,iC,ic2,nC)]=1;
                            if (stateprm[jp]==1 && fabs(p1-1)<OTOLERANCE && fabs(mp->tpm21-1)<OTOLERANCE) stateprm[iK(iC,iC,ic2,nC)]=2;
                        }
                        if (statepr[jp] && fabs(p2-1)<OTOLERANCE && stateprm[jp]!=1 && mp->pti[iC*nC+iC]<OTOLERANCE)
                        {   //only if technology certainly improves leave a trace of who had the move on this level
                            stateprm[iK(iC,ic1,iC,nC)]=stateprm[jp];
                        }
                        if (statepr[jp] && fabs(p2-1)<OTOLERANCE && stateprm[jp]!=1 && iC==0)
                        {
                            if (stateprm[jp]==2 && fabs(p2-1)<OTOLERANCE && fabs(mp->tpm12-1)<OTOLERANCE) stateprm[iK(iC,ic1,iC,nC)]=1;
                            if (stateprm[jp]==2 && fabs(p2-1)<OTOLERANCE && fabs(mp->tpm22-1)<OTOLERANCE) stateprm[iK(iC,ic1,iC,nC)]=2;
                        }
                    }
                    else
                    {   //simultaneous move
                        statepr[iK(iC,iC,ic2,nC)]=statepr[iK(iC,iC,ic2,nC)] || (statepr[jp] && p1>0);
                        statepr[iK(iC,ic1,iC,nC)]=statepr[iK(iC,ic1,iC,nC)] || (statepr[jp] && p2>0);
                        statepr[iK(iC,iC,iC,nC)] =statepr[iK(iC,iC,iC,nC)]  || (statepr[jp] && p1>0 && p2>0);
                        //if (statepr[jp]) printf("interior iC=%d ic1=%d ic2=%d jp=%d (p1=%1.5f>0)=%d (p2=%1.5f>0)=%d\n",iC,ic1,ic2,jp,p1,(p1>0),p2,(p2>0));
                    }
                }
            }
            //SECOND look at edges
            //Assume: technology always has a positive probability of improvement (so order of moves in current layer can not be traced)
            //The only exception is iC==0
            for (ic1=iC+1;ic1<nC;ic1++)
            {
                jp=iK(iC,ic1,iC,nC); //current point
                p1=g[iC].p1[id(nC,ic1,iC,iC,eqstring[jp])];
                if (mp[0].alternate) statepr[iK(iC,iC,iC,nC)]=statepr[iK(iC,iC,iC,nC)] || (statepr[jp] && (p1>0) && stateprm[jp]!=2);
                else                 statepr[iK(iC,iC,iC,nC)]=statepr[iK(iC,iC,iC,nC)] || (statepr[jp] && (p1>0));
                if (statepr[jp] && iC==0 && mp[0].alternate)
                {
                    if (stateprm[jp]==1 && fabs(p1-1)<OTOLERANCE && fabs(mp->tpm11-1)<OTOLERANCE) stateprm[iK(iC,iC,iC,nC)]=1;
                    if (stateprm[jp]==1 && fabs(p1-1)<OTOLERANCE && fabs(mp->tpm21-1)<OTOLERANCE) stateprm[iK(iC,iC,iC,nC)]=2;
                }
            }
            for (ic2=iC+1;ic2<nC;ic2++)
            {
                jp=iK(iC,iC,ic2,nC); //current point
                p2=g[iC].p2[id(nC,iC,ic2,iC,eqstring[jp])];
                if (mp[0].alternate) statepr[iK(iC,iC,iC,nC)]=statepr[iK(iC,iC,iC,nC)] || (statepr[jp] && (p2>0) && stateprm[jp]!=1);
                else                 statepr[iK(iC,iC,iC,nC)]=statepr[iK(iC,iC,iC,nC)] || (statepr[jp] && (p2>0));
                if (statepr[jp] && iC==0 && mp[0].alternate)
                {
                    if (stateprm[jp]==2 && fabs(p2-1)<OTOLERANCE && fabs(mp->tpm12-1)<OTOLERANCE) stateprm[iK(iC,iC,iC,nC)]=1;
                    if (stateprm[jp]==2 && fabs(p2-1)<OTOLERANCE && fabs(mp->tpm22-1)<OTOLERANCE) stateprm[iK(iC,iC,iC,nC)]=2;
                }
            }
            //THIRD look at corner : except that there are no state-changing investment decisions in the corner
        }
        // printeqb(statepr,neqstr,iC,NULL); //print again after the investment probs are accounted for
        // if (mp[0].alternate) printeqb(stateprm,neqstr,90+iC,NULL); //print again after the investment probs are accounted for
        
        //mark lower layers using technology progress pti
        if (iC>0)
        {
            for (ic1=iC;ic1<nC;ic1++)
            {
                for (ic2=iC;ic2<nC;ic2++)
                {
                    if (statepr[iK(iC,ic1,ic2,nC)])
                    {   
                        jp=iK(iC,ic1,ic2,nC); //current point
                        p1=g[iC].p1[id(nC,ic1,ic2,iC,eqstring[jp])];
                        p2=g[iC].p2[id(nC,ic1,ic2,iC,eqstring[jp])];
                        if (mp[0].alternate)
                        {   //alternating move
                            if (   ic1==iC || ic2==iC  //corner and edges
                                || ((p1<1 || stateprm[iK(iC,ic1,ic2,nC)]==2) && (p2<1 || stateprm[iK(iC,ic1,ic2,nC)]==1)) //interior point where game could stay for another period
                                )
                            {   //fill points below if positive prob to fall through to there
                                for (jC=1;jC<=iC;jC++)
                                {
                                    jp=(iC-jC)*nC+iC;
                                    statepr[iK(iC-jC,ic1,ic2,nC)]=statepr[iK(iC-jC,ic1,ic2,nC)] || (mp->pti[jp]>0);
                                    //in case of deterministic technology, check for deterministic alternation
                                    if (fabs(mp->pti[jp]-1)<OTOLERANCE && stateprm[iK(iC,ic1,ic2,nC)]>0)
                                    {
                                        if (stateprm[iK(iC,ic1,ic2,nC)]==1 && fabs(mp->tpm11-1)<OTOLERANCE) stateprm[iK(iC-jC,ic1,ic2,nC)]=1;
                                        if (stateprm[iK(iC,ic1,ic2,nC)]==1 && fabs(mp->tpm21-1)<OTOLERANCE) stateprm[iK(iC-jC,ic1,ic2,nC)]=2;
                                        if (stateprm[iK(iC,ic1,ic2,nC)]==2 && fabs(mp->tpm12-1)<OTOLERANCE) stateprm[iK(iC-jC,ic1,ic2,nC)]=1;
                                        if (stateprm[iK(iC,ic1,ic2,nC)]==2 && fabs(mp->tpm22-1)<OTOLERANCE) stateprm[iK(iC-jC,ic1,ic2,nC)]=2;
                                    }
                                }
                                //finish up with current layer (jC=0)
                                jp=iC*nC+iC;
                                if ((ic1==iC || ic2==iC) && iC<nC-1)
                                {
                                    statepr[iK(iC,ic1,ic2,nC)]=(mp->pti[jp]>OTOLERANCE); //rewrite mark on edges if technology certainly improves
                                    stateprm[iK(iC,ic1,ic2,nC)]=0; //because technology may stay or to delete temporary mark which is not processed to lower layers
                                } 
                            }
                        }
                        else
                        {   //simultaneous move
                            //in the interior if p==1 the point should not be translated downwards
                            if (ic1==iC || ic2==iC || (p1<1-OTOLERANCE && p2<1-OTOLERANCE))
                            {   //fill points below if positive prob to be there
                                for (jC=0;jC<=iC;jC++)
                                {
                                    jp=(iC-jC)*nC+iC;
                                    if (jC==0 && iC<nC-1 && (ic1==iC || ic2==iC)) statepr[iK(iC-jC,ic1,ic2,nC)]=(mp->pti[jp]>0); //rewrite mark on edges if technology certainly improves
                                    else statepr[iK(iC-jC,ic1,ic2,nC)]=statepr[iK(iC-jC,ic1,ic2,nC)] || (mp->pti[jp]>0);

                                }
                            }
                        } 
                    }
                }
            }
        }
    }//next iC
    //statepr[iK(iC,ic1,ic2,nC)] is computed! Point (iC,ic1,ic2) may applear on equilibrium path IFF statepr[iK(iC,ic1,ic2,nC)] == true
    //stateprm[iK(iC,ic1,ic2,nC)] for alternating move game is also computed! For simultaneous move ==0

    //check feasibility indicators
    // printeqb(statepr,neqstr,1,NULL);
    // printeqb(stateprm,neqstr,99,NULL);

    pure=1; //pure strategies in ALL FEASIBLE nodes
    symmetry=1; //symmetry in ALL nodes
    leapfrog=0; //leapfrogging (positive prob of high cost to invest) in at least some FEASIBLE nodes
    for (iC=0;iC<nC;iC++)
    {
        for (ic1=iC;ic1<nC;ic1++)
        {
            for (ic2=iC;ic2<nC;ic2++)
            {
                ieqb=eqstring[iK (iC,ic1,ic2,nC)];//selected stage equilibrium
                
                //pure strategy when investment probs for both are 0 or 1 when FEASIBLE and RIGHT TO MOVE
                if (  statepr[iK(iC,ic1,ic2,nC)] && 
                      ((fabs(g[iC].p1[IdOUT]-0)>OTOLERANCE && stateprm[iK(iC,ic1,ic2,nC)]!=2 &&  fabs(g[iC].p1[IdOUT]-1)>OTOLERANCE) || 
                       (fabs(g[iC].p2[IdOUT]-0)>OTOLERANCE && stateprm[iK(iC,ic1,ic2,nC)]!=1 && fabs(g[iC].p2[IdOUT]-1)>OTOLERANCE))
                    ) pure=0;
                //symmetric when value functions are symmetric
                //NB Need to look at selected equilibrium for reveresed ic1/ic2
                if (mp[0].alternate==1)
                {   //alternate move
                    if (   fabs(g[iC].v10[id(nC,ic1,ic2,iC,ieqb)]-g[iC].v20[id(nC,ic2,ic1,iC, eqstring[iK(iC,ic2,ic1,nC)] )])>OTOLERANCE ||
                           fabs(g[iC].v11[id(nC,ic1,ic2,iC,ieqb)]-g[iC].v21[id(nC,ic2,ic1,iC, eqstring[iK(iC,ic2,ic1,nC)] )])>OTOLERANCE ||
                           fabs(g[iC].v20[id(nC,ic1,ic2,iC,ieqb)]-g[iC].v10[id(nC,ic2,ic1,iC, eqstring[iK(iC,ic2,ic1,nC)] )])>OTOLERANCE ||
                           fabs(g[iC].v21[id(nC,ic1,ic2,iC,ieqb)]-g[iC].v11[id(nC,ic2,ic1,iC, eqstring[iK(iC,ic2,ic1,nC)] )])>OTOLERANCE ||
                           fabs(g[iC].x10[id(nC,ic1,ic2,iC,ieqb)]-g[iC].x20[id(nC,ic2,ic1,iC, eqstring[iK(iC,ic2,ic1,nC)] )])>OTOLERANCE ||
                           fabs(g[iC].x11[id(nC,ic1,ic2,iC,ieqb)]-g[iC].x21[id(nC,ic2,ic1,iC, eqstring[iK(iC,ic2,ic1,nC)] )])>OTOLERANCE ||
                           fabs(g[iC].x20[id(nC,ic1,ic2,iC,ieqb)]-g[iC].x10[id(nC,ic2,ic1,iC, eqstring[iK(iC,ic2,ic1,nC)] )])>OTOLERANCE ||
                           fabs(g[iC].x21[id(nC,ic1,ic2,iC,ieqb)]-g[iC].x11[id(nC,ic2,ic1,iC, eqstring[iK(iC,ic2,ic1,nC)] )])>OTOLERANCE
                        ) symmetry=0;
                }
                else
                {   //simultanious move
                    if (   fabs(g[iC].v10[id(nC,ic1,ic2,iC,ieqb)]-g[iC].v20[id(nC,ic2,ic1,iC, eqstring[iK(iC,ic2,ic1,nC)] )])>OTOLERANCE ||
                           fabs(g[iC].v11[id(nC,ic1,ic2,iC,ieqb)]-g[iC].v21[id(nC,ic2,ic1,iC, eqstring[iK(iC,ic2,ic1,nC)] )])>OTOLERANCE ||
                           fabs(g[iC].v20[id(nC,ic1,ic2,iC,ieqb)]-g[iC].v10[id(nC,ic2,ic1,iC, eqstring[iK(iC,ic2,ic1,nC)] )])>OTOLERANCE ||
                           fabs(g[iC].v21[id(nC,ic1,ic2,iC,ieqb)]-g[iC].v11[id(nC,ic2,ic1,iC, eqstring[iK(iC,ic2,ic1,nC)] )])>OTOLERANCE
                        ) symmetry=0;
                }
                //leapfrogging when probability to invest is positive when FEASIBLE and RIGHT TO MOVE
                if (  statepr[iK(iC,ic1,ic2,nC)] &&
                      ((ic1>ic2 && g[iC].p1[IdOUT]>0 && stateprm[iK(iC,ic1,ic2,nC)]!=2) ||
                       (ic1<ic2 && g[iC].p2[IdOUT]>0 && stateprm[iK(iC,ic1,ic2,nC)]!=1))
                    )
                    {
                        leapfrog=1;
                    }

                //underinvestment
                if (mp->monsoluiton!=NULL)
                {   //only for when gmon is passed in last (optional) argument
                    //CONDITION: monopoly invests AND neither of firms invests with prob==1
                    if (    (mp->monsoluiton[iK(iC,MIN(ic1,ic2),nC-1,nC)]>0 ) &&   
                            (fabs(g[iC].p1[id(nC,ic1,ic2,iC,ieqb)]-1.0)>OTOLERANCE) && 
                            (fabs(g[iC].p2[id(nC,ic1,ic2,iC,ieqb)]-1.0)>OTOLERANCE)
                        )
                    {
                        underinvestment++; //sum of points where the condition
                        // printf("test %d %f\n",(fabs(g[iC].p1[id(nC,ic1,ic2,iC,ieqb)]-1.0)<OTOLERANCE),g[iC].p1[id(nC,ic1,ic2,iC,ieqb)]-1.0+OTOLERANCE);
                        // printf("underinvestment mon=%d p1=%f p2=%f\n",mp->monsoluiton[iK(iC,MIN(ic1,ic2),nC-1,nC)],g[iC].p1[id(nC,ic1,ic2,iC,ieqb)],g[iC].p2[id(nC,ic1,ic2,iC,ieqb)]);
                    }   
                }
                if (pure==0 && symmetry==0 && leapfrog==1) break;
            }
            if (pure==0 && symmetry==0 && leapfrog==1) break;
        }
        if (pure==0 && symmetry==0 && leapfrog==1) break;
    }

    // check for underinvestment (compare to monopoly solution)
    underinvestment=0; //indicator for underinvestment compared to monopoly
    if (mp->monsoluiton!=NULL) //only fill out the map in case 1 equilibrium is to be computed
    {
        for (iC=0;iC<nC;iC++)
        {
            for (ic1=iC;ic1<nC;ic1++)
            {
                for (ic2=iC;ic2<nC;ic2++)
                {
                    ieqb=eqstring[iK (iC,ic1,ic2,nC)];//selected stage equilibrium
                    j=iK(iC,MIN(ic1,ic2),nC-1,nC);
                    //j=iK(iC,ic1,ic2,nC);
                    //Fill out the map for current ESR, THIRD solumn => 2*
                    mp->monvsown[j] =  (mp->monvsown[j]>0) || 
                         (fabs(g[iC].p1[id(nC,ic1,ic2,iC,ieqb)]-1.0)<OTOLERANCE) || 
                         (fabs(g[iC].p2[id(nC,ic1,ic2,iC,ieqb)]-1.0)<OTOLERANCE); /////////////////CONDITION FOR THE INDICATOR in 4th row
                    //printf("iC=%d (%d,%d) min=%d iK=%d p1=%f p2=%f ind=%f\n",iC,ic1,ic2,MIN(ic1,ic2),j,g[iC].p1[id(nC,ic1,ic2,iC,ieqb)],g[iC].p2[id(nC,ic1,ic2,iC,ieqb)],mp->monvsown[j]);
                }
            }
            //compare to monopoly
            for (ic1=iC;ic1<nC;ic1++)
            {
                for (ic2=iC;ic2<nC;ic2++)
                {
                    j=iK(iC,MIN(ic1,ic2),nC-1,nC);
                    if (mp->monsoluiton[j]>0 && fabs(mp->monvsown[j])<OTOLERANCE) 
                    {
                        underinvestment=1;
                        break;
                    }
                }
                if (underinvestment==1) break;
            }
            if (mp->esrmax>1) break; //only break if we are not filling out the full map
        }
    }
#endif

    //PART 2 Save the output variables
    outindex0=outindex;//save value of outindex before adding
    //mandatory variables to be ALWAYS saved in this order
    output[outindex++]=(double) lexindex;  // lexicographic index
    output[outindex++]=1;                  // initial number of duplicates of this point (only one is outputed)
    

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   VARIABLES TO OUTPUT
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//IMPORTANT : NUMBER OF ADDITIONAL VALUES FOR OUTPUT MUST BE CORRECT, otherwise CRASH !!!!!!
//            If these are changed, change labels in graph.EqbstrPlot()!
#define NUMOUTVALUES 7
    output[outindex++]=g[nC-1].v10[0]; // value : not investing : firm1
    if (mp[0].alternate) output[outindex++]=g[nC-1].x20[0]; // value : not investing : firm2, turn of firm 1, only one stage equilibrium
    else                 output[outindex++]=g[nC-1].v20[0]; // value : not investing : firm2
    output[outindex++]=pure;           // type of eqb path
    output[outindex++]=symmetry;       // type of eqb path
    output[outindex++]=leapfrog;       // type of eqb path
    // efficiency score
    if (mp[0].alternate) output[outindex++]=g[nC-1].ecm1[0]; // value: firm1 making the first move (consistent with above)
    else                 output[outindex++]=g[nC-1].ecm1[0]; // value: simultanious move game stored in ecm1
    output[outindex++]=underinvestment;  // type of eqb path

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //PART 3 Check for uniquiness if not first output
    if (outindex0>0) 
    {   //number of elements already outputed is
        length0=outindex0/(NUMOUTVALUES+2);
        //search for the new element among already saved
        bsres=bsearch(output+outindex0,output,length0,sizeof(double)*(NUMOUTVALUES+2),&compout);
        //printf("searching for\n<%f %f %f %f %f>\namong\n",*(output+outindex0),*(output+outindex0+1),*(output+outindex0+2),*(output+outindex0+3),*(output+outindex0+4));
        //for (i=0;i<length0;i++) 
        //printf("[%f %f %f %f %f]\n",*(output+i*(NUMOUTVALUES+2)),*(output+i*(NUMOUTVALUES+2)+1),*(output+i*(NUMOUTVALUES+2)+2),*(output+i*(NUMOUTVALUES+2)+3),*(output+i*(NUMOUTVALUES+2)+4));
        if (bsres)
        {   //bsres!=NULL ==> such values already exist!
            //update the counter of duplicates 
            bsres[1]++;//second element
            //forget about new outputstring
            outindex = outindex0;//reset value for the next call 
            //cleanup
            for (i=0;i<NUMOUTVALUES+2;output[outindex0+i++]=mxGetNaN());
            //printf("duplicate\n");
        }
        else
        {   //bsres==NULL ==> the value is new
            //resort the output
            //printf("saving and re-sorting\n");
            qsort(output,length0+1,sizeof(double)*(NUMOUTVALUES+2),&compout);
        }
    }
    //else printf("saved first\n");
    
    //print eqb string
    #if PRINTeqbstr==1
// CHOOSE WHAT EQBSTRINGS TO OUTPUT HERE    
//    if (floor(g[nC-1].v10[0])==36 && floor(g[nC-1].v20[0])==1) 
//    if (floor(g[nC-1].v10[0])==36) 
//    if (g[nC-1].ecm1[0] > 60.5198 - 1e-4 )
// if (g[nC-1].ecm1[0] > 86.4192 - 1e-4)
        printeqb(eqstring,neqstr,lexindex,statepr);
    #endif


#ifdef OUTSTATISTICS
    free(statepr);
    free(stateprm);
#endif

    return 0;
}

static int compout(const void *a,const void *b) 
{   //comparing function to be used in output routine
    //compares blocks of NUMOUTVALUES+2 doubles using by 2-3-4-.. columns using tolerance
    int i;
    for (i=2;i<NUMOUTVALUES+2;i++)
    {
        if      (*((double *)a+i) > *((double *)b+i) +OTOLERANCE) return 1;
        else if (*((double *)b+i) > *((double *)a+i) +OTOLERANCE) return -1;
        //here a[i] and b[i] are equal, go to next digit 
        //printf("%f == %f\n",*((double *)a+i),*((double *)b+i));    
    }
    return 0; //equality
}

void lexistr(int *eqstring,int nstr,size_t index)
{   //fills out the space eqstr with lexicographical string of given index
    //uses MAXEQB as module 
    size_t p1=MAXEQB,p2=1,i;
    for (i=0;i<nstr;i++)
    {
        eqstring[i]=(index%p1)/p2;
        p1*=MAXEQB;
        p2*=MAXEQB;
    }
}

size_t lexindex(double *eqstring,int nstr)
{   //returns the index of passed eqstring
    //uses MAXEQB as module
    size_t res=0, m=1; 
    int i;
    for (i=0;i<nstr;i++)
    {
        res+=(int)eqstring[i]*m;
        m*=MAXEQB;
    }
    return res;
}

size_t lexindex1(int *eqstring,int nstr)
{   //returns the index of passed eqstring
    //uses MAXEQB as module
    size_t res=0,m=1; 
    int i;
    for (i=0;i<nstr;i++)
    {
        res+=eqstring[i]*m;
        m*=MAXEQB;
    }
    return res;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int prerun(MPstruct mp, Gamestruct* g, Bnestruct bne, Brstruct* br, double* OutESRtable,int* eqstring1,size_t *indx) {
/*--------------------------------------------------------------------------------------
 This function runs the model up to a given start string to check for eligibility
 on each integer. In the end eqstr contains the starting esr string and g structure contains
 the results of a run with this esr
 Output: 1 if start eqstring infeasible, 0 if feasible
         eqstring1 contains last feasible ESR from the cycle
**--------------------------------------------------------------------------------------*/
    int j,nC;
    int neqstr;
    int iCz=0,ic1z=0,ic2z=0;
    int *eqstring0;
    //init dimentions
    nC=mp.nC;
    neqstr=(nC*(nC+1)*(2*nC+1)/6); //number of digits in eqstring
    //allocate the space for the temp equilibruim selection strings
    eqstring0=(int *) calloc(neqstr,sizeof(int));
    //initialize eqsting0
    if (mp.esrstartpt) for (j=0;j<neqstr;eqstring0[j]=(int)mp.esrstartpt[neqstr-1-j++]);//firstr eqbstring was passed
    else lexistr(eqstring0,neqstr,mp.esrstart);//the first eqbstring to run 
    if (PRINTeqbloop>0) printf("Start eqbstring:\n");
    if (PRINTeqbloop>0) printeqb(eqstring0,neqstr,lexindex1(eqstring0,neqstr),NULL);
    if (PRINTeqbloop>0) printf("Checking for eligibility starting from higher digits:\n");
    //initialize eqsting1 to zeros
    for (j=0;j<neqstr;eqstring1[j++]=0);
    //initialize g structure by solving the model for zeros
    if (mp.alternate==1) leap_am(mp, g, bne, br, eqstring1, ic1z, ic2z, iCz);
    else leap_sm(mp, g, bne, br, eqstring1, ic1z, ic2z, iCz);
    //run the cycle to add digits one by one and check feasibility on each step
    for (j=neqstr;j>=0;j--)
    {   //run from higher digits <=> lower layer of the game
        if (eqstring0[j]!=eqstring1[j])
        {   //skip identical digits
            // 1 check feasibility for this digit
            iKinv (j,nC,&iCz,&ic1z,&ic2z);
//            if (eqstring0[j]>=g[iCz].neqb[id(nC,ic1z,ic2z,iCz,0)]) break;//not feasible
            if (eqstring0[j]>=g[iCz].neqb[id(nC,ic1z,ic2z,iCz,0)])
            {
            printf("Infeasibility at digit %d which is %d, g.neqb(iC=%d,ic1=%d,ic2=%d)=%d\n",j,eqstring0[j],iCz,ic1z,ic2z,g[iCz].neqb[id(nC,ic1z,ic2z,iCz,0)]);

             mexWarnMsgTxt("Initial ESR string passed to the solver is infeasible, returning.");
             break;//not feasible
            }
            // 2 replace next digit
            eqstring1[j]=eqstring0[j];
            if (PRINTeqbloop>0) printeqb(eqstring1,neqstr,lexindex1(eqstring1,neqstr),NULL);
            // 3 recalculate g structure
            // initialize the g structure from icz level up
            
            if (mp.alternate==1) leap_am(mp, g, bne, br, eqstring1, ic1z, ic2z, iCz);
            else leap_sm(mp, g, bne, br, eqstring1, ic1z, ic2z, iCz);
        }
    }
    free(eqstring0);//free
    indx[0]=lexindex1(eqstring1,neqstr);
    if (j<0) return 0; //all feasible!
    else return 1; //error: infeasible start
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void run(MPstruct mp, Gamestruct* g, Bnestruct bne, Brstruct* br, double* OutESRtable) {
/*--------------------------------------------------------------------------------------
 MAIN "do the work" function started from mexFunction
**--------------------------------------------------------------------------------------*/
    size_t index,i,ii;
    int j,jj,nC;
    int neqstr;
    int iCz=0,ic1z=0,ic2z=0;
    int *eqstring;
    int outerr=0;
    
    nC=mp.nC;
    outputiter(0,NULL,NULL,NULL,NULL); //initializing call to output routine
    
    //announce the type of game
    if (PRINTeqbloop>0) {
        if (mp.alternate==1) printf("Solving alternating moves game\n");
        else printf("Solving simutanious moves game\n");
    }
    
    if (mp.esr!=ESRstring)
    {   //run only ONCE using classic ESR
        //NB: ESRstring representation of the rule will be outputed in forth output variable
        if (mp.alternate==1) leap_am(mp, g, bne, br, NULL, ic1z, ic2z, iCz);
    	    else leap_sm(mp, g, bne, br, NULL, ic1z, ic2z, iCz);
    }
    else
    {   //run the CYCLE OVER eqstrings
        neqstr=(nC*(nC+1)*(2*nC+1)/6); //number of digits in eqstring
        //allocate the space for the equilibruim selection strings
        eqstring=(int *) calloc(neqstr,sizeof(int));
        for (j=0;j<neqstr;eqstring[j++]=0);
        //initialize counters
        i=0;//the index of eqstring in lexicographical order (0 to MAXEQB^neqstr)
        ii=0;//number of feasible eqilibrium paths examined
        j=0;//the digit last changed in the lexicographical steps
        //restrict the number of runs here!
        while (true) //cycle over eqstring lexicographical generation
        {
            // 1 step: do the work
            //solve
            if (mp.esrstart>0 && i==0)
            {   //run special cycle if starting eqbstirng is not zero
                if (prerun(mp, g, bne, br, OutESRtable,eqstring,&i)==0)
                {
                    ii++;
                    outerr=outputiter(i,eqstring,&mp,g,OutESRtable);
                }
                else return; //exit if infeasible 
                //else mexWarnMsgTxt("Initial ESR string passed to the solver is infeasible.");
                if (mp.esrmax==1) break;//special case for 1 output
            }
            else
            {
                ii++;//counts the runs of solver
                if (mp.alternate==1) leap_am(mp, g, bne, br, eqstring, ic1z, ic2z, iCz);
                else leap_sm(mp, g, bne, br, eqstring, ic1z, ic2z, iCz);
                //output
                outerr=outputiter(i,eqstring,&mp,g,OutESRtable);
            }
            //output the eqstirng if asked for
            if (PRINTeqbloop>1) printeqb(eqstring,neqstr,i,NULL);
    
            //2 step: make one lexicographical step
            eqstring[0]++;
            i++;//update i
            for (j=0;j<neqstr-1;j++) //run mod arithmetics
            {
                if (eqstring[j]>=MAXEQB)
                {
                    eqstring[j]=0;
                    eqstring[j+1]++;
                }
                else break;
            } //on exit digit j is the only one that changed up
            
            //3 step: check feasibility condition for the updated digit
            //        and jump over the infeasible strings
            //calculate iC,ic1,ic2
            iKinv (j,nC,&iCz,&ic1z,&ic2z);
            //also check for the last eqstring
            while (eqstring[neqstr-1]<MAXEQB && eqstring[j]>=g[iCz].neqb[id(nC,ic1z,ic2z,iCz,0)])
            {
                if (PRINTeqbloop>1)
                    printf("Jumping over %1.0f eqstrings because of infeasibility\n",(MAXEQB-1)*pow(MAXEQB,(double)j));
                i+=(MAXEQB-1)*(int)pow(MAXEQB,(double)j);//update i
                //concerned about OVERFLOW in i
                //if (i>OVERFLOW) break;
                //make the jump
                if (j>=neqstr-1) break;
                eqstring[j+1]++;//increase next digit by 1
                for (jj=j+1;jj<neqstr-1;jj++) //run mod arithmetics
                {
                    if (eqstring[jj]>=MAXEQB)
                    {
                        eqstring[jj]=0;
                        eqstring[jj+1]++;
                    }
                    else break;
                } //on exit digit jj is the only one that changed up
                while (j>=0) eqstring[j--]=0;//set j and all digits below j to 0
                //update the last changed index
                j=jj;
                iKinv (j,nC,&iCz,&ic1z,&ic2z);
            }
            //5 step: check if no more eqstrings
            if (eqstring[neqstr-1]>=MAXEQB || (j==neqstr-1 && eqstring[neqstr-1]>=g[iCz].neqb[id(nC,ic1z,ic2z,iCz,0)])) break;
        }//end of while over eqstrings
        if (PRINTeqbloop>0 || outerr!=0) printf("Checked %u equilibrium selection rules and found %d equilibria(um)\n",i,ii);
        //free
        free(eqstring);
    }
}

mxArray * monopoly_solution_map(MPstruct *mp,const mxArray * gmon)
{   //this function constructs the maps of monopoly investments
    mxArray *out, *layer;
    double *dout, *dlayer;
    int iC,ic1,ic2,j,i;
    int nC,niK,nm;

    nC=mp->nC;
    niK=nC*(nC+1)*(2*nC+1)/6;
    //create matrix for output to Matlab
    out=mxCreateDoubleMatrix(niK,4,mxREAL);
    dout=(double*) mxGetPr(out);

    for (j=0;j<niK*4;dout[j++]=mxGetNaN());
    for (iC=0;iC<nC;iC++)
    {
        for (ic1=iC;ic1<nC;ic1++)
        {
            for (ic2=iC;ic2<nC;ic2++)
            {
                j=iK(iC,MIN(ic1,ic2),nC-1,nC); //map ic1,ic2 into min()
                //j=iK(iC,ic1,ic2,nC); //full index for testing
                // STRUCTURE OF THE ADDITIONAL output matrix
                dout[0*niK+j]=iC;           // 1 iC
                dout[1*niK+j]=MIN(ic1,ic2); // 2 MIN
                dout[2*niK+j]=0;//init 0!   // 3 Monopoly investments
                dout[3*niK+j]=0;//init 0!   // 4 current ESR investrments
                //make available to internal C function
                mp->monsoluiton=dout+(nC*(nC+1)*(2*nC+1)/6)*2; //monopoly investment indicators will be here
                mp->monvsown=   dout+(nC*(nC+1)*(2*nC+1)/6)*3; //own investment indicators should be put here                
            }
        }
        //PARSE G structure of monopoly
        // ASSUMPTIONS:
        // Gcol colums 13 - selected equilibrium
        //              1 - ic1
        //              2 - ic2
        //             11 - p1 investment probability of firm1:=monopoly
        layer=mxGetField(gmon,iC,"solution");
        dlayer=(double*) mxGetPr(layer);
        if (dlayer==NULL) mexErrMsgTxt("error parsing the passed game structure!");
        nm=mxGetM(layer);
        for (i=0;i<nm;i++)
        {
            if (!mxIsNaN(dlayer[(13-1)*nm+i]) && fabs(dlayer[(13-1)*nm+i]-1.0)<OTOLERANCE)
            {   //chosen equilibrium
                ic1=(int) dlayer[(1-1)*nm+i];
                ic2=(int) dlayer[(2-1)*nm+i];
                j=iK(iC,MIN(ic1,ic2),nC-1,nC);
                //j=iK(iC,ic1,ic2,nC); //full index for testing
//printf("iC=%d (%d,%d) min=%d iK=%d pr=%f\n",iC,ic1,ic2,MIN(ic1,ic2),j,dlayer[(11-1)*nm+i]);
                dout[2*niK+j]= dout[2*niK+j] || (fabs(dlayer[(11-1)*nm+i]-1.0)<OTOLERANCE);  //////////////////// CONDITION FOR MONOPOLY MAP INDICATOR
            }
        }
    }

    //return pointer to created Matlab matrix
    return out;
}




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const char *keys[] = {"solution","c","ec"};
    const int nkeys = 3;
    const char *brkeys[] = {"br"};
    const int nbrkeys = 1;
    Gamestruct *g;
    Bnestruct bne;
    Brstruct *br;
    MPstruct mp;
    double *OutESRtable;
    int j, iC, nC;
    
    if (nrhs != 3 && nrhs != 4) mexErrMsgTxt("Error in leapfrog.c: wrong number of inputs!");
    if (nlhs != 4 && nlhs != 5) mexErrMsgTxt("Error in leapfrog.c: wrong number of outputs!");
    if ((mxGetM(prhs[0])!=1) && (mxGetN(prhs[1])!=9)) mexErrMsgTxt("Wrong parameter vector, params!");
    if ((mxGetM(prhs[1])!=1) && (mxGetN(prhs[1])!=7)) mexErrMsgTxt("Wrong parameter vector, modparams!");
    if ((mxGetM(prhs[2])!=1) && (mxGetN(prhs[2])!=5)) mexErrMsgTxt("Wrong switches vector, sw!");
    
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
    mp.pti     =            mxGetPr(mxGetField(prhs[0],0,"pti"));
    
    //   model parameters
    //   (see setup.m for description mp)
    mp.dt       =   mxGetScalar(mxGetField(prhs[1],0,"dt"));
    mp.df       =   mxGetScalar(mxGetField(prhs[1],0,"df"));
    mp.k1       =   mxGetScalar(mxGetField(prhs[1],0,"k1"));
    mp.k2       =   mxGetScalar(mxGetField(prhs[1],0,"k2"));
    mp.sigma    =   mxGetScalar(mxGetField(prhs[1],0,"sigma"));
    mp.c_og     =   mxGetScalar(mxGetField(prhs[1],0,"c_og"));
    mp.eta      =   mxGetScalar(mxGetField(prhs[1],0,"eta"));
    mp.c_tr     =   mxGetScalar(mxGetField(prhs[1],0,"c_tr"));
    mp.tpm     =    mxGetPr(mxGetField(prhs[1],0,"tpm"));
    mp.tpm11=mp.tpm[0]; mp.tpm21=mp.tpm[2];   
    mp.tpm12=mp.tpm[1]; mp.tpm22=mp.tpm[3];

    
    //   switches controling the
    //   (see setup.m for description mp)
    mp.alternate=   mxIsLogicalScalarTrue(mxGetField(prhs[2],0,"alternate"));
    mp.analytical=   mxIsLogicalScalarTrue(mxGetField(prhs[2],0,"analytical"));
    mp.esr      =   (int)mxGetScalar(mxGetField(prhs[2],0,"esr"));
    mp.esrstartpt = NULL; //initialize
    if (mxGetM(mxGetField(prhs[2],0,"esrstart"))==1 && mxGetN(mxGetField(prhs[2],0,"esrstart"))==1) //scalar=index of lexistring
    {
        mp.esrstart =   (size_t)mxGetScalar(mxGetField(prhs[2],0,"esrstart"));
        if (mp.esrstart > pow(MAXEQB,nC*(nC+1)*(2*nC+1)/6)-1) mexErrMsgTxt("The passed index for the initial eqbstring is too large!");
    }
    else if (mxGetN(mxGetField(prhs[2],0,"esrstart"))!=(nC*(nC+1)*(2*nC+1)/6) && mxGetM(mxGetField(prhs[2],0,"esrstart"))!=1)
        mexErrMsgTxt("sw.esrstart contains a vector of wrong dimensions for the starting ESR string!");
    else
    {
        mp.esrstartpt = (double*) mxGetPr(mxGetField(prhs[2],0,"esrstart")); //pointer to passed eqbstr to start with
        mp.esrstart = lexindex(mp.esrstartpt,(nC*(nC+1)*(2*nC+1)/6));
    }
    mp.esrmax   =   (size_t)mxGetScalar(mxGetField(prhs[2],0,"esrmax"));
    mp.esrmaxout = mp.esrmax*(NUMOUTVALUES+2); //max number of output lines (for convenience)

    
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
        mxSetField(plhs[2], iC, keys[0], mxCreateDoubleMatrix((nC-iC)*(nC-iC)*MAXEQB, 20, mxREAL));
        mxSetField(plhs[2], iC, keys[1], mxCreateDoubleScalar(mxGetNaN()));
        mxSetField(plhs[2], iC, keys[2], mxCreateDoubleMatrix((nC-iC)*(nC-iC)*MAXEQB, 10, mxREAL));
        
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
        //Additional values for alternating move game
        //the values when it's other firm's turn to invest
        g[iC].x10=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*15;
        g[iC].x11=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*16;
        g[iC].x20=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*17;
        g[iC].x21=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*18;
        g[iC].neqb=g[iC].ic1+((nC-iC)*(nC-iC)*MAXEQB)*19;
		//Expected costs
        g[iC].ecv10=(double *) mxGetData(mxGetField(plhs[2], iC, keys[2])); //key= ec
        g[iC].ecv11=g[iC].ecv10+((nC-iC)*(nC-iC)*MAXEQB)*1;
        g[iC].ecv20=g[iC].ecv10+((nC-iC)*(nC-iC)*MAXEQB)*2;
        g[iC].ecv21=g[iC].ecv10+((nC-iC)*(nC-iC)*MAXEQB)*3;
        g[iC].ecx10=g[iC].ecv10+((nC-iC)*(nC-iC)*MAXEQB)*4;
        g[iC].ecx11=g[iC].ecv10+((nC-iC)*(nC-iC)*MAXEQB)*5;
        g[iC].ecx20=g[iC].ecv10+((nC-iC)*(nC-iC)*MAXEQB)*6;
        g[iC].ecx21=g[iC].ecv10+((nC-iC)*(nC-iC)*MAXEQB)*7;
        g[iC].ecm1 =g[iC].ecv10+((nC-iC)*(nC-iC)*MAXEQB)*8;
        g[iC].ecm2 =g[iC].ecv10+((nC-iC)*(nC-iC)*MAXEQB)*9;

		//Initialize to NaNs
        for (j=0;j<(nC-iC)*(nC-iC)*MAXEQB*(20-1);g[iC].ic1[j++]=mxGetNaN()); //leave last column without NaNs
        for (j=0;j<(nC-iC)*(nC-iC)*MAXEQB*10;g[iC].ecv10[j++]=mxGetNaN());
        g[iC].c=(double *) mxGetData(mxGetField(plhs[2], iC, keys[1])); //key= c
   
    }

    /********************************************************/
    // OUTPUT: table with esr cycle results
    /********************************************************/

    if (mp.esr!=ESRstring)
    {   //if ESR is not eqstring
        //output eqstring implied by the ESR
        plhs[3]=mxCreateDoubleMatrix(1,(nC*(nC+1)*(2*nC+1)/6),mxREAL);
        mp.out4=(double*) mxGetData(plhs[3]);
    }
    else
    {   //if ESR is eqstring, output the eq table
        plhs[3]=mxCreateDoubleMatrix(NUMOUTVALUES+2,mp.esrmax,mxREAL);
        OutESRtable=(double*) mxGetData(plhs[3]);
        for (j=0;j<mp.esrmax*(NUMOUTVALUES+2);OutESRtable[j++]=mxGetNaN());
    }


    /********************************************************/
    // OUTPUT: additional output: comparison to monopoly
    /********************************************************/
    //Optional input parameter: gamestructure of solution of monopoly case (for checking current solutions against monopoly)
    if (nrhs==4)
    {
        if (nlhs<5) mexErrMsgTxt("Not enough output arguments, expecting 5 arguments for monopoly comparison case");
        else
        {
            if (mxIsStruct(prhs[3]) && mxGetNumberOfElements(prhs[3])==nC)
            {   //consutruc map of statespace and mark where monopoly makes investment decisions
                plhs[4]=monopoly_solution_map(&mp,prhs[3]);
            }
            else mexErrMsgTxt("Unexpected optional argument, expecting game structure of the same dimention which is interpreted as monopoly solution");
        }
    }
    else
    {
        mp.monsoluiton=NULL; //monopoly investment indicators would be here if monopoly game strucutre was passed
        mp.monvsown=NULL;

    }    

    run(mp, g, bne, br, OutESRtable);
    //call the main part of the model
    
    //cleanup
    free(mp.cgrid);
    free(g);
    free(br);
}







