//         Fedor Iskhakov, University Technology Sidney
//         John Rust, University of Maryland
//         Bertel Schjerning, University of Copenhagen
//         February 2011

// Method for looking for equilibria: EQBMETHOD_GRID, EQBMETHOD_VF
#define EQBMETHOD_GRID

//Possible ESRs: ESR_HighCostInvests ESR_MixedStr ESR_FirstInvests ESR_HighCostInvestsMix
#define ESR_HighCostInvestsMix
// #define ESR_HighCostInvests



double f_kf(int iC, MPstruct *mp) {
    /*----------------------------------------------------------------------------------------------------------------
     *          Cost of technological improvement
     **----------------------------------------------------------------------------------------------------------------*/
    return mp[0].k1/(1+mp[0].k2*mp[0].cgrid[iC]);
}

double f_pti(int iC, MPstruct *mp) {
    /*----------------------------------------------------------------------------------------------------------------
     *          Probability of technological improvement
     **----------------------------------------------------------------------------------------------------------------*/
    if (iC==0) return 0.0; // no way to improve from c=0
    return mp[0].c_tr*mp[0].cgrid[iC]/(1.0+mp[0].c_tr*mp[0].cgrid[iC]);
}

int f_neqb(int nC, int ic1, int ic2, int iC, Gamestruct g[]) {
    /*----------------------------------------------------------------------------------------------------------------
     *          Number of equilibria
     **----------------------------------------------------------------------------------------------------------------*/
    int i, j;
    //find number of eqb for these c1,c2,c
    i=id(nC,ic1,ic2,iC,0);
    for (j=0;j<MAXEQB;j++) {
        if (mxIsNaN(g[iC].ieqb[i+j])) break;
    }
//printf("c1c2c=%d-%d-%d Number of equilibria : %d\n",ic1,ic2,iC,j);
    return j;
}

int f_ers(int iC, int ic1, int ic2, int neqb) {
    /*----------------------------------------------------------------------------------------------------------------
     *          Equilibrium selection rules
     **----------------------------------------------------------------------------------------------------------------*/
// ESR_HighCostToInvest
//         equilibrium selection rule for the  leapfrogging game
//         This function specifies an index for the selected equilibrium
//         in every state of the game. For states where there is only one
//         equilibrium, the function will return the value 0, and if there
//         are n equilibria, it will return an integer index for the selected
//         equilibria between  0 and neqb-1. Usually neqb is 3 and the equilibria are
//         ordered with neqb=1 corresponding to the "no investment for firm 1"
//         equilibrium (probability equal or close to 0 that firm 1 invests),
//         neqb=1 is an interior mixed strategy equilibrium, and neqb=3 is the "investment
//         by firm 1" equilibrium (probability equal or close to 1 that firm 1 invests)
    
    // only one equilibrium at the edges FOR ANY ESR
    if ((ic1 == iC) || (ic2 == iC)) return 0;
    
    // if equilibrium is unique - select it
    if (neqb == 1) return 0;
    
#ifdef ESR_HighCostInvests
//ESR1 - similar to John
// The firm with the highest cost gets to invest and in case of a tie, firm 1 gets to invest
if (ic1 < ic2) return 0;          // firm 2 gets to invest
if (ic1 > ic2) return neqb-1;     // firm 1 gets to invest
if (ic1 == ic2) return neqb-1;    // firm 1 gets to invest in case of tie
#endif
#ifdef ESR_HighCostInvestsMix
//ESR1 - similar to John
// The firm with the highest cost gets to invest and in case of a tie, firm 1 gets to invest
if (ic1 < ic2) return 0;          // firm 2 gets to invest
if (ic1 > ic2) return neqb-1;     // firm 1 gets to invest
//     if (ic1 > ic2) return f_neqb(nC, ic1, ic2, iC, g)-1;     // firm 1 gets to invest
if (ic1 == ic2) { // pick mixed equilibrium when possible
    return 1;
//         if (f_neqb(nC, ic1, ic2, iC, g)==3) return 1;
//         if (f_neqb(nC, ic1, ic2, iC, g)==5) return 2;
//         if (f_neqb(nC, ic1, ic2, iC, g)==7) return 3;
    return 0;
}
#endif
#ifdef ESR_MixedStr
//This eqb rule always picks mixed equilibrium when possible
if (neqb==3) return 1;
if (neqb==5) return 2;
if (neqb==7) return 3;
return 0;
#endif
#ifdef ESR_FirstInvests
//Highest prob of investment of firm 1 (last equilibrium)
return neqb-1;
#endif
//rise an error if nothing was returned by now
mexErrMsgTxt("Equilibrium selection rules not chosen or no equilibrium selected for given values of state variables!");

}

void s_bne(MPstruct *mp, Bnestruct *bne) {
    /*--------------------------------------------------------------------------------------
     *        solves the "end game" equilibria when costs using the state of the
     *        art production technology have reached the lowest possible level, which is normalized
     *        to zero and are treated as an absorbing state
     **--------------------------------------------------------------------------------------*/
    
    
    // NOTE: s_bne does not appropriately take outside good into account when sigam>0
    
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
    // BN-equaliprium prices, profits and market shares
    /********************************************************/
    for (ic1=0;ic1<nC;ic1++) { //ic1
        for (ic2=0;ic2<=ic1;ic2++) {  //ic2
            c1=mp[0].cgrid[ic1];    // very much used => for easy referencing
            c2=mp[0].cgrid[ic2];
            
            // Pure Bertrand solution
            if (mp[0].og==0) { //no outside good
                f_logit(nF, mp[0].sigma, s, (-1)*c1, (-1)*c2); //put logit probs into pi
                p[0]=MAX(c1, c2);
            }
            else { //outside good
                f_logit3(nF+mp[0].og, mp[0].sigma, s, (-1)*c1, (-1)*c2, (-1)*mp[0].c_og); //put logit probs into pi
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
                if (mp[0].og>0) mexErrMsgTxt("Don't know how to handle sigma>0 and par.og>0!");
                cp=2*mp[0].ctol;
                p[2]=p[0];
                /*find prices that satisfy F.O.C. of prifit maximization
                 * for both firms => solve 2x2 system of unlinear equations
                 * by Newton method:
                 * i=1..3: sigma-(1-s[i])(p[i]-c[i])=0, where s[i] is logit prob of full x
                 * Bertrand prices are used as startting values */
                
                for (it=0;((it<mp[0].maxit) && (cp>=mp[0].ctol));it++) { //iterations
                    //f_logit3(nF+mp[0].og, mp[0].sigma, s, (-1)*p[0], (-1)*p[1], (-1)*p[2]); //put logit probs into pi
                    f_logit(nF, mp[0].sigma, s, (-1)*p[0], (-1)*p[1]); //put logit probs into pi
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
                //no probabilities here
            }
            else mexErrMsgTxt("Don't know how to handle sigma<0!");
            pf[0]=(p[0]-c1)*s[0];
            pf[1]=(p[1]-c2)*s[1];
            
            // Save output
            i=ic1*nC+ic2;
            j=ic2*nC+ic1;
            
            bne[0].ic1[i]=ic1;
            bne[0].ic2[i]=ic2;
            bne[0].c1[i]=mp[0].cgrid[ic1];
            bne[0].c2[i]=mp[0].cgrid[ic2];
            
            bne[0].p1[i]=p[0];
            bne[0].p2[i]=p[1];
            bne[0].pf1[i]=pf[0];
            bne[0].pf2[i]=pf[1];
            bne[0].s1[i]=s[0];
            bne[0].s2[i]=s[1];
            
            // by symetry
            bne[0].ic1[j]=ic2;
            bne[0].ic2[j]=ic1;
            bne[0].c1[j]=mp[0].cgrid[ic2];
            bne[0].c2[j]=mp[0].cgrid[ic1];
            
            bne[0].p2[j]=bne[0].p1[i];
            bne[0].p1[j]=bne[0].p2[i];
            bne[0].pf2[j]=bne[0].pf1[i];
            bne[0].pf1[j]=bne[0].pf2[i];
            bne[0].s2[j]=bne[0].s1[i];
            bne[0].s1[j]=bne[0].s2[i];
        }
    }
}


