/*        John Rust, University of Maryland
 *        Bertel Schjerning, University of Copenhagen
 *        Fedor Iskhakov, University Technology Sidney
 *        November 2011
 *
 *  EXPERIMENTAL new findeqb algorithm
 *  Main objective: speed and independence of small adjustments (as compared to Bertel's original algorithm)
*/
#include "stdio.h"
#include <stdlib.h>
#include "math.h"
#include "matrix.h"
#include "mex.h"
#include <omp.h>
#define MAX(A,B) ((A)>(B)?(A):(B))
#define MIN(A,B) ((A)>(B)?(B):(A))
#define PI 3.14159265358979323846264338327950288419716939

#ifndef DBL_MAX
#define DBL_MAX 1e300
#endif

#define MAXEQB 500

typedef struct {
    int maxit; 
    double ctol;
    double bounds[2];
} MPstruct; //structure of parameters

double testfun(double x,double* bounds)
{   //returns the value of the function we test
    //unless x is outside of given region
    mxArray *in, *out;
    //return inf outside of true bounds 
    if (x<bounds[0]) return -DBL_MAX;
    if (x>bounds[1]) return  DBL_MAX;
    //else return func value
    in = mxCreateDoubleScalar(x);
    mexCallMATLAB(1,&out,1,&in,"testfun");
    return mxGetScalar(out);
}

double sa (double start,double* lb,double* ub,MPstruct* mp)
{   //successive approximations up to exit from the interval [lb,ub]
    //returns stable fixed point or DBL_MAX/-DBL_MAX
    int it;
    double tmp,tmp1;
    for (it=0;it<mp[0].maxit;it++)
    {
        tmp=testfun(start,mp[0].bounds);
printf(" sa: it=%d bounds=[%1.17f,%1.17f] testfun(%1.17f)=%1.17f\n",it,lb[0],ub[0],start,tmp);        
        if (tmp<lb[0]) return -DBL_MAX;
        if (tmp>ub[0]) return DBL_MAX;
        //use usual criterion as first step
        if (fabs(tmp-start)<mp[0].ctol)
        {   //stop if tmp and tmp+ctol are on different sides of fixed point
printf(" sa: stopping1: tmp-start=%1.5e approaching from %s\n",tmp-start,(tmp-start>0?"below":"above"));
            if (tmp-start>0)
            {   //approaching from below
                tmp1=testfun(tmp+mp[0].ctol,mp[0].bounds);
                if (tmp1-tmp-mp[0].ctol<0)
                {   //full stop
printf(" sa: stopping2: tmp-start=%1.5e tmp1-tmp-ctol=%1.5e\n",tmp-start,tmp1-tmp-mp[0].ctol);
                    lb[0]=tmp+mp[0].ctol;//update lb so that lb is to the right of fx point
printf(" sa: result: tmp=%1.17f, lb=%1.17f, lb-tmp=%1.5e\n",(tmp+tmp1)/2,lb[0],lb[0]-(tmp+tmp1)/2);
                    return (tmp+tmp1)/2;
                }
            }
            else
            {   //approaching from above
                tmp1=testfun(tmp-mp[0].ctol,mp[0].bounds);
                if (tmp1-tmp+mp[0].ctol>0)
                {   //full stop
printf(" sa: stopping2: tmp-start=%1.5e tmp1-tmp+ctol=%1.5e\n",tmp-start,tmp1-tmp+mp[0].ctol);
                    lb[0]=tmp;//update lb so that lb is to the right of fx point
printf(" sa: result: tmp=%1.17f, lb=%1.17f, lb-tmp=%1.5e\n",(tmp+tmp1)/2,lb[0],lb[0]-(tmp+tmp1)/2);
                    return (tmp+tmp1)/2;
                }
            }
        }
        start=tmp;
    }
    mexWarnMsgTxt("Max number of iterations exceeded in sa!");
}

double br0(MPstruct *mp,double* domain)
{   //finds fixed point on the domain by successive bisections
        double u,l,m,uf,lf,mf;
        int it=0;
        l=domain[0];
        u=domain[1];
        lf=testfun(l,mp[0].bounds)-l;
        uf=testfun(u,mp[0].bounds)-u;
        //check if there is fixed point at all at this segment
        if ((lf<0 && uf<0) || (lf>0 && uf>0))  mexErrMsgTxt("Incorrect call of bracketing algorithm!");
        //search for fixed point
        while (fabs(u-l)>mp[0].ctol && it<=mp[0].maxit)
        {
            m=(u+l)/2;
            mf=testfun(m,mp[0].bounds)-m;
            if ((mf>0 && uf>0) || (mf<0 && uf<0)) u=(u+l)/2;
            else l=(u+l)/2;
            it++;
        }
        if (it>=mp[0].maxit) mexWarnMsgTxt("Max number of iterations exceeded in br0!");
        return (u+l)/2;
}


void findfxpts(double *out,MPstruct* mp,double* domain)
{   //this function finds fixed points on testfun
    int it,outi=0; //indeces for iteration and output
    double tmp, l, u; //tmp value and upper bound 
    double lb=domain[0], ub=domain[1]; //lower and upper bound

    //1. check I that eqbs exist
    //we could also update upper bound in case of stable eqb,
    //but this is only ok for monotonically increasing testfun
printf("check I\n");
    tmp=sa(domain[1],&lb,&ub,mp);//lb is updated!
    if (tmp==-DBL_MAX) return;//no equilibria

    //2. check II that eqbs exist
printf("check II\n");
    lb=domain[0];
    tmp=sa(domain[0],&lb,domain+1,mp);//lb is updated!
    if (tmp==DBL_MAX) return;//no equilibria

    //3. continue with the main cycle
    //tmp and lb already initialized!
    while (tmp<DBL_MAX && outi<MAXEQB)
    {
        if (tmp==-DBL_MAX)
        {   //next fixed point is unstable ==>
printf("tmp=-DBL_MAX last_eqb=%1.16f lb=%1.16f ub=%f Bracketing\n",out[outi-1],lb,ub);
            //do the bracketing to find point where tmp!=-DBL_MAX
            //which will be an unstable fixed point
            it=0;
            u=ub;
            while (fabs(l-u)>mp[0].ctol && it<mp[0].maxit)
            {
                it++;
//                if (sa((lb+u)/2,&lb,&u,mp)>-DBL_MAX) u=(lb+u)/2;
//                else lb=(lb+u)/2; //permanently update lower bound
                l=lb;//different l to be updated in sa
                l=lb+mp[0].ctol;
                if (sa((lb+u)/2,&l,&u,mp)>-DBL_MAX)
                {
                    u=(lb+u)/2;
printf("BR: it=%d func>-DBL_MAX ==> lb=%1.16f u=%f \n",it,lb,u);                
                }
                else
                {
                    lb=(lb+u)/2; //permanently update lower bound
printf("BR: it=%d  func=-DBL_MAX ==> lb=%1.16f u=%f \n",it,lb,u);
                }
            }
            if (it>=mp[0].maxit) mexWarnMsgTxt("Max number of iterations exceeded in unstable eqb search!");
            //unless eld of domain save unstable equilibrium
            if (fabs(u-ub)>mp[0].ctol) out[outi++]=(lb+u)/2;
printf("unstable eqb=%f found after %d iterations, outi=%d\n",out[outi-1],it,outi);
            //step over unstable equilibrium
            lb=u;
        }
        else
        {   //next fixed point is stable
printf("tmp=%f  lb=%f ub=%f (lb>tmp)=%d\n",tmp,lb,ub,(lb>tmp));
            //save and step over
            out[outi++]=tmp;
printf("stable eqb=%f, outi=%d\n",out[outi-1],outi);
            //lb is already update in sa (to step a little over tmp)
        }
        //check ctol in the end of the cycle
        if (fabs(lb-ub)<mp[0].ctol) break;
        //look at next segment
printf("end of cycle");
        tmp=sa(lb,&lb,&ub,mp);
printf("tmp=%f\n",tmp);

    }
    if (outi>=MAXEQB) mexWarnMsgTxt("Max number of equilibria exceeded!");
}

void getdomain(double* domain,MPstruct* mp) 
{   //this function finds "unknown" domain with given tolerance
    //Input: domain = initial search area
    //Output: domain = exact domain
        double u, l, tmp;
        int it;

        l=domain[0];
        u=domain[1];
        //search for max
        if (testfun(u,mp[0].bounds)==DBL_MAX)
        {   //search only if there is space above, otherwise keep the upper bound
            it=0;
            while (fabs(u-l)>mp[0].ctol && it<=mp[0].maxit)
            {
                if (testfun((u+l)/2,mp[0].bounds)<DBL_MAX) l=(u+l)/2;
                else u=(u+l)/2;
                it++;
            }
            if (it>=mp[0].maxit) mexWarnMsgTxt("Max number of iterations exceeded in domain (max)!");
            domain[1]=l;
        }

        l=domain[0];
        u=domain[1];
        //search for min
        if (testfun(l,mp[0].bounds)==-DBL_MAX)
        {   //search only if there is space below, otherwise keep the lower bound
            it=0;
            while (fabs(u-l)>mp[0].ctol && it<=mp[0].maxit)
            {
                if (testfun((u+l)/2,mp[0].bounds)>-DBL_MAX) u=(u+l)/2;
                else l=(u+l)/2;
                it++;
            }
            if (it>=mp[0].maxit) mexWarnMsgTxt("Max number of iterations exceeded in domain (min)!");
            domain[0]=u;
        }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{   //Gateway: should be called with initial bounds (between 0 and 1), returns a list of equilibria
    double bounds[2],*out;
    int i;
    MPstruct mp;
    double domain[2];

    if (nrhs != 2) mexErrMsgTxt("Error: wrong number of inputs!");
    if (nlhs != 1) mexErrMsgTxt("Error: wrong number of outputs!");
    if (mxGetNumberOfElements(prhs[0]) != 2) mexErrMsgTxt("Error: input must be [lowerbound, upperbound]!");
    //read in input
    mp.bounds[0]=*((double*)mxGetPr(prhs[0]));
    mp.bounds[1]=*((double*)mxGetPr(prhs[0])+1);
    mp.ctol=mxGetScalar(prhs[1]); //tolerance
    mp.maxit=1000;   //max nr of iterations
    //create output
    plhs[0]=mxCreateDoubleMatrix(1,MAXEQB,mxREAL);
    out=(double*) mxGetPr(plhs[0]);
    for (i=0;i<MAXEQB;out[i++]=mxGetNaN());

    //1 find domain
    domain[0]=0.0;//initial search region
    domain[1]=1.0;
    getdomain(domain,&mp);
    printf("true bounds = [%2.15f,%2.15f]\n",mp.bounds[0],mp.bounds[1]);
    printf("     domain = [%2.15f,%2.15f]\n",domain[0],domain[1]);
    printf("max abs domain error = %1.5e\n",MAX(fabs(domain[0]-mp.bounds[0]),fabs(domain[1]-mp.bounds[1])));
    printf("                ctol = %1.5e\n",mp.ctol);

    //2 find equilibria
    findfxpts(out,&mp,domain);

}

