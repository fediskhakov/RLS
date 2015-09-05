// void s_domainbr2(Brcstruct *brc,MPstruct *mp,double* domain);
// double s_findeqb_bracketing(Brcstruct *brc,MPstruct *mp,double* domain);
// int sa(double* x1, double x0, double l, double r,Brcstruct* brc,MPstruct *mp);
// void findeqb_sa(Brcstruct* brc,MPstruct *mp,MVstruct *mv,double* domain);
// double f_invbr(double p_i, Brcstruct* brc);

static int VERBOSE;

void analytical_am(MPstruct* mp, Brc_amstruct* brc, Gamestruct* gc, int ic1, int ic2, int iC, int nC) {
    /*--------------------------------------------------------------------------------------
     * This function calculates all equilibria for given c1 c2 c
     * and puts them in proper output structures gc and brc
     *--------------------------------------------------------------------------------------*/
    int i, j, iF, idz;
    int ieqb, neqb;
    int save;
    double a0[2], a1[2], a2[2], a3[2], b0[2], b1[2], c1[2], c2[2], f[2], pstar[4][2], A[2], B[2], v0_0[2], v0_1[2], pbr[2];
    double D0[2], D1[2], D2[2];
    int nPstar[2];
    
    // ********************************************
    // Coefficients
    // ********************************************
    //precompute for convenience
    
    a0[0]=brc[0].pf[0]+brc[0].h11;
    a1[0]=brc[0].pf[0]+brc[0].h12;
    a2[0]=brc[0].betnpc*mp[0].tpm21*brc[0].x1[0];
    a3[0]=brc[0].betnpc*mp[0].tpm22*brc[0].x1[0];

    a0[1]=brc[0].pf[1]+brc[0].h22;
    a1[1]=brc[0].pf[1]+brc[0].h21;
    a2[1]=brc[0].betnpc*mp[0].tpm12*brc[0].x1[1];
    a3[1]=brc[0].betnpc*mp[0].tpm11*brc[0].x1[1];
    
    b0[0]=brc[0].betnpc*mp[0].tpm11;
    b1[0]=brc[0].betnpc*mp[0].tpm12;
    b0[1]=brc[0].betnpc*mp[0].tpm22;
    b1[1]=brc[0].betnpc*mp[0].tpm21;

    
    c1[0]=brc[0].betnpc*mp[0].tpm21;
    c2[0]=brc[0].betnpc*mp[0].tpm22;
    c1[1]=brc[0].betnpc*mp[0].tpm12;
    c2[1]=brc[0].betnpc*mp[0].tpm11;
    

    for (iF=0;iF<2;iF++) {
        D0[iF]=c1[iF]*a1[iF]
                +(1-c2[iF])*a0[iF]
                +brc[0].v1[iF]*(c2[iF]+b0[iF]+c1[iF]*b1[iF]-c2[iF]*b0[iF]-1);
        D1[iF]=a2[iF]
                +c1[iF]*(a3[iF]-a1[iF])
                +c2[iF]*(a0[iF]-a2[iF])
                +(c2[iF]*b0[iF]-c2[iF]-c1[iF]*b1[iF])*brc[0].v1[iF];
        D2[iF]= -(c1[iF]*a3[iF]-c2[iF]*a2[iF]);
        D2[iF]= 0.0;    // True for alternating move game

    }
           
    // ********************************************
    // Switch by eta
    // ********************************************
    if (mp[0].eta==0.0) {   //Simple case when best response function have constant coefficients
        // Equilibrium candidates, NOTE iF index is the opponent here !!!!!!
        for (iF=0;iF<2;iF++) {
            if (fabs(D2[iF])<1e-10) {
                pstar[0][iF]=0.0;
                pstar[1][iF]=-D0[iF]/D1[iF]; // root
                pstar[2][iF]=1;
                pstar[3][iF]=-99; // XXX
                
                nPstar[iF]=3;
            }
            else {
                // Equilibrium candidates, NOTE iF index is the opponent herD   !!!!!!
                pstar[0][iF]=0.0;
                pstar[1][iF]=-(1/(2*D2[iF]))*(D1[iF]+sqrt(D1[iF]*D1[iF]-4*D0[iF]*D2[iF]));//positive root
                pstar[2][iF]=-(1/(2*D2[iF]))*(D1[iF]-sqrt(D1[iF]*D1[iF]-4*D0[iF]*D2[iF]));//negative root
                pstar[3][iF]=1;
                nPstar[iF]=4;
            }
        }
        
        //check all combinations of roots (and corners)
        ieqb=-1;
        save=0;
        for (j=0;j<nPstar[1];j++) {
            for (i=0;i<nPstar[0];i++) {
                //loop twice over pstar
                // check four corners of the reaction function box for pure strategy equilibria
                if (((i==0) || (i==(nPstar[0]-1))) && ((j==0)  || (j==(nPstar[1]-1)))) { // corner
                    //best responce for firm 1:
                    pbr[0]=((D0[0]+D1[0]*pstar[i][0])<0);
                    //best responce for firm 2:
                    pbr[1]=((D0[1]+D1[1]*pstar[j][1])<0);
                    //check if probs really match
                    if (fabs(pbr[0]-pstar[j][1])<mp[0].ctol && fabs(pbr[1]-pstar[i][0])<mp[0].ctol){
                        ieqb=ieqb+1;
                        save=1;
                    }
                }
                if (((i==1) || (i==(nPstar[0]-2))) && ((j==1) || j==(nPstar[1]-2))) { // mixed strategy
                    
                    if (((pstar[i][0]>=0) && (pstar[i][0]<=1)) && ((pstar[j][1]>=0) && (pstar[j][1]<=1))) {
                        ieqb=ieqb+1;
                        save=2;
                    }
                }
                if (ieqb>=MAXEQB) mexErrMsgTxt("Error in analytical_am.c, eta==0: number of equilibria exceeded MAXEQB");
                //save new equilibrium (note IsOUT depends on ieqb)
                if (save>0) {
                    
                    idz=IdOUT;
                    gc[0].ieqb[idz]=ieqb;          // already populated
                    gc[0].eqbtype[idz]=(save==1?0:1); //0 = pure, 1 = mixed
                    gc[0].v10[idz]=brc[0].v0[0];
                    gc[0].v11[idz]=brc[0].v1[0];
                    gc[0].v21[idz]=brc[0].v1[1];
                    gc[0].p1[idz]=pstar[j][1];
                    gc[0].p2[idz]=pstar[i][0];
                    gc[0].x11[idz]=brc[0].x1[0];
                    gc[0].x21[idz]=brc[0].x1[1];
                    
                    //simplify
                    f[0]=c1[0]*(1-gc[0].p2[idz])/(1-c2[0]*(1-gc[0].p2[idz]));
                    B[0]=b0[0]+b1[0]*f[0];
                    A[0]=a0[0]+a1[0]*f[0]+a2[0]*gc[0].p2[idz]+a3[0]*gc[0].p2[idz]*f[0];
                    
                    f[1]=c1[1]*(1-gc[0].p1[idz])/(1-c2[1]*(1-gc[0].p1[idz]));
                    B[1]=b0[1]+b1[1]*f[1];
                    A[1]=a0[1]+a1[1]*f[1]+a2[1]*gc[0].p1[idz]+a3[1]*gc[0].p1[idz]*f[1];
                    
                    if (save==1) { //pure strategy
                        //simplify again
                        v0_1[0]=A[0]+B[0]*gc[0].v11[idz];
                        v0_1[1]=A[1]+B[1]*gc[0].v21[idz];
                        for (iF=0;iF<2;iF++) {
                            v0_0[iF]=A[iF]/(1-B[iF]);
                        }
                        gc[0].v10[idz]=v0_1[0]*pbr[0]+v0_0[0]*(1-pbr[0]);
                        gc[0].v20[idz]=v0_1[1]*pbr[1]+v0_0[1]*(1-pbr[1]);
                    }
                    else { //mixed strategy
                        gc[0].v10[idz]=A[0]+B[0]*gc[0].v11[idz];
                        gc[0].v20[idz]=A[1]+B[1]*gc[0].v21[idz];
                    }
                    gc[0].x10[idz]=(brc[0].pf[0]+brc[0].betnpc*mp[0].tpm12*f_lnsum(2, mp[0].eta, gc[0].v11[idz],gc[0].v10[idz])
                        +brc[0].betnpc*mp[0].tpm22*gc[0].p2[idz]*brc[0].x1[0]+brc[0].h12)
                        /(1-brc[0].betnpc*mp[0].tpm22*(1-gc[0].p2[idz]));
                    gc[0].x20[idz]=(brc[0].pf[1]+brc[0].betnpc*mp[0].tpm21*f_lnsum(2, mp[0].eta, gc[0].v21[idz],gc[0].v20[idz])
                        +brc[0].betnpc*mp[0].tpm11*gc[0].p1[idz]*brc[0].x1[1]+brc[0].h21)
                        /(1-brc[0].betnpc*mp[0].tpm11*(1-gc[0].p1[idz]));
                    
                    
                    
                    //reset the save flag
                    save=0;
                    

                }
            }
        }
        //the number of found equilibria
        neqb=ieqb+1;
        gc[0].neqb[id(nC,ic1,ic2,iC,0)]=neqb;
    }
    else {  //eta>0
        mexErrMsgTxt("Error in analytical_am.c, Analytical solution not implemented for eta>0: ");
    }
}


