//         Fedor Iskhakov, University Technology Sidney
//         John Rust, University of Maryland
//         Bertel Schjerning, University of Copenhagen
//         February 2011
#include <stdarg.h>




double rmax(int n, double *arr) {//returns the max of the array
/*----------------------------------------------------------------------------------------------------------------
 * Recursive max function: Computes the max of over all n elements in the array arr
-----------------------------------------------------------------------------------------------------------------*/
    double maybemax;
    if (n==1) return *arr;
    maybemax = rmax(n-1, arr+1);
    if (*arr>=maybemax) return *arr;
    else return maybemax;
}

double f_lnsum(int nr,double eta, ...) {
/*----------------------------------------------------------------------------------------------------------------
** f_lnsum: 
 * Computes social surpus function for the logit (LogSum)
 *
 * inputs:  j, number of alternatives 
 * outputs: 
**----------------------------------------------------------------------------------------------------------------*/
    va_list ap;
    double mxm, v, sum;
    int i;
    
    va_start(ap, eta);
    mxm=va_arg(ap, double);
    for (i=1;i<nr;i++){
        v=va_arg(ap, double);
        if (v>mxm) mxm=v;
    }
    
    if (eta<=0.0){
        //max
        return mxm;
    }
    else {
        //real log-sum
        va_start(ap, eta);
        sum=0;
        va_start(ap, eta);
        for (i=0;i<nr;i++){
            sum=sum+exp((va_arg(ap, double)-mxm)/eta);
        }

        va_end(ap);
        return mxm+eta*log(sum);
    }
}

double f_logit(int nr, double sigma, double *pr, ...) {
//----------------------------------------------------------------------------------------------------------------
    /*Calculates logit probability for the first option of nr `utility` levels with scaling parameter sigma
      Also if *pr (size nr) is passed, it is filled with probabilities for all options
      Variable argument list are the nr "utility" levels
      Sigma=0 is correctly processed as max 
     */
    va_list ap;
    double ave=0, sum=0, frst, mxm, v;
    int i, j=0;
    
    va_start(ap, pr);
    mxm=va_arg(ap, double);
    frst=mxm;
    for (i=1;i<nr;i++){
        v=va_arg(ap, double);
        if (v>mxm) mxm=v;
    }

    if (sigma<=0.0){
        //max ==> prob=1 or 0
        //count number of args==max
        va_start(ap, pr);
        for (i=0;i<nr;i++, j+=(mxm==va_arg(ap, double))); //j is number of max

        //calculate probs
        va_start(ap, pr);
        if (pr!=NULL) for(i=0;i<nr;pr[i++]=(mxm==va_arg(ap, double)?1/(double)j:0.0));
        va_end(ap);
        return ((mxm==frst)?(1/(double)j):0.0);
    }
    else {
        //logit probs
        sum=0;
        va_start(ap, pr);
        for (i=0;i<nr;i++, sum+=exp((va_arg(ap, double)-mxm)/sigma));
        va_start(ap, pr);
        if (pr!=NULL) for(i=0;i<nr;pr[i++]=(exp((va_arg(ap, double)-mxm)/sigma))/sum);
        va_end(ap);
        return (exp((frst-mxm)/sigma))/sum;
    }
}


double f_solvevf(double a, double b, double c, double eta, MPstruct *mp){
    int it;
    double v,v1,cp=1.0;
//  Purpose: Solve scalar value function equations on the form v=a+b*logsum(v,c) w.r.t. v using Newtons method
//  inputs:    a,b,c:  constants in value function
//             eta:    parameter in logsum function
//             v0:     starting value
//  output:    v:      solution v=a+b*logsum(v,c)

    //good starting point - from asymptotic approximation of logsum when v goes to infinity   
    v=a*eta/(1-b*eta);
    for (it=0;((it<mp[0].maxit) && (cp>=mp[0].tolerance));it++){
        //Newton step
        v1 = v-(v-a-b*f_lnsum(2, eta, v, c))/(1-b*f_logit(2, eta, NULL, v, c));
//        printf("it=%d v=%g v1=%g ", it, v, v1);
        cp=fabs(v1-v);
//        printf("tol=%g\n", tol);
        v=v1;
    }
    if (it>=mp[0].maxit) mexErrMsgTxt("Failed to converge in v=a+b*logsum(v,c) solution (see function.c/f_solvef)!");
    return (v1);
}


