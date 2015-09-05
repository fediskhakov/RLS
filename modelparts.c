//         Fedor Iskhakov, University Technology Sidney
//         John Rust, University of Maryland
//         Bertel Schjerning, University of Copenhagen
//         February 2011



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
    if (mp[0].c_tr<0) { // Deterministic technological improvement each period until c=0
        return 1.0;
    }
    else { // STOCHASTIC IMPORVEMENT
        return mp[0].dt*mp[0].c_tr*mp[0].cgrid[iC]/(1.0+mp[0].dt*mp[0].c_tr*mp[0].cgrid[iC]);
        
        
    }
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
            bne[0].pf1[i]=pf[0]*mp[0].dt;
            bne[0].pf2[i]=pf[1]*mp[0].dt;
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


