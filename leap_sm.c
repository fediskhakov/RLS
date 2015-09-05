#define EULER 0.577215664901532


void leap_sm(MPstruct mp, Gamestruct* g, Bnestruct bne, Brstruct* br, int* eqstring, int ic1z, int ic2z, int iCz) {
    /********************************************************/
// SECTION 1: DECLARATIONS
    /********************************************************/
    // Params and Modparams
    int nC; //Main parameter
    MVstruct mv;
    
    //double v, v1, r, p1;
    int i, j, iC, jp, jC, ic1, ic2, ieqb, neqb;
    int i1, i2, iccc;
    //double q; //temps for c>0
    //double v, v1, r, p1;
    //int i, j, iF, iI, iC, ic1, ic2, ieqb, it, ip, neqb;
    double v, q;
    double* hc1;
    double* hc2;
    Brc_amstruct brc;
    double EP10, EP11, EP20, EP21, logsum1, logsum2, p1, p2, f1;//constants for symbolic solutions in cost recursions
    double AN1,AI1,B,AN2,AI2;
#ifdef PRINTnrstagempe    
    int neqstr; //only for output of ESR for a given number of stage equilibria
#endif
    
    nC        =  mp.nC;               // very much used => for easy referencing
    
    /********************************************************/
// SECTION 2: MEMORY ALLOCATION
    /********************************************************/
    mv.h1      = (double *) calloc(nC*nC, sizeof(double));     // temp storage of formula parts, precomputing, see c>0
    mv.h2      = (double *) calloc(nC*nC, sizeof(double));     // temp storage of formula parts, precomputing, see c>0
    hc1        = (double *) calloc(nC*nC, sizeof(double));     // temp storage of formula parts, precomputing, see c>0
    hc2        = (double *) calloc(nC*nC, sizeof(double));     // temp storage of formula parts, precomputing, see c>0
    mv.eqbinfo = (double *) calloc(5*MAXEQBINFO, sizeof(double));     // storage for found equilibria in findeqb
    
    
    
    /********************************************************/
// SECTION 3: LOOP OVER iC
    /********************************************************/
    for (iC=iCz;iC<nC;iC++) {
        /********************************************************/
        // Initialize for iC
        
        /********************************************************/
        g[iC].c[0]=mp.cgrid[iC]; // save value of c in output
        
        mv.pc=f_pti(iC,&mp);    // probability of a technological improvement. Initially improvements are assumed to be
        // 'incremental', that is the jump, if it occurs, is to c'=cgrid(i-1) with probability pti
        // for iC=0 pti=0!
        mv.kv=f_kf(iC, &mp);     // Investment costs
        mv.logsumK=f_lnsum(2, mp.eta, 0, -mv.kv);        // Update mv.kv using investments costs, K(c)
        
        // precompute h(c1,c2,c) first element of H(c1,c2,c);
        if (iC==0) for (i=0;i<nC*nC;i++) {mv.h1[i]=0;mv.h2[i]=0; hc1[i]=0;hc2[i]=0;}
        else {
//        omp_set_num_threads(NUM_THREADS);
//        #pragma omp parallel shared(mv, g, mp,iC,  nC) private(ic1, ic2, neqb, i, j)
            {
//               printf("There are %d threads\n",omp_get_num_threads());
//        #pragma omp for
                
                for (ic1=iC;ic1<nC;ic1++) {
                    for (ic2=iC;ic2<nC;ic2++) {

                        j=(ic1-iC)*(nC-iC)+(ic2-iC);    // index used for storage in h
                        mv.h1[j]=0;
                        mv.h2[j]=0;
                        hc1[j]=0;
                        hc2[j]=0;
                        for (jC=1;jC<=iC;jC++) {
                            jp=(iC-jC)*nC+iC;
                            if (mp.pti[jp]>0) {
                                neqb=(int)g[iC-jC].neqb[id(nC,ic1,ic2,iC-jC,0)];
                                i=id(nC, ic1, ic2, iC-jC, f_ers(&mp, iC-jC, ic1, ic2, neqb,eqstring));    // index for selected equilibrium
                                mv.h1[j]=mv.h1[j]+mp.pti[jp]*f_lnsum(nF, mp.eta, g[iC-jC].v10[i], g[iC-jC].v11[i]);
                                mv.h2[j]=mv.h2[j]+mp.pti[jp]*f_lnsum(nF, mp.eta, g[iC-jC].v20[i], g[iC-jC].v21[i]);
                                hc1[j]=hc1[j]+mp.pti[jp]*(g[iC-jC].p1[i]*g[iC-jC].ecv11[i]+(1-g[iC-jC].p1[i])*g[iC-jC].ecv10[i]);  
                                hc2[j]=hc2[j]+mp.pti[jp]*(g[iC-jC].p2[i]*g[iC-jC].ecv21[i]+(1-g[iC-jC].p2[i])*g[iC-jC].ecv20[i]);
//                                printf("f_ers(&mp, iC-jC=%d, ic1=%d, ic2=%d, neqb=%d,eqstring)=%d\n",iC-jC, ic1, ic2, neqb, f_ers(&mp, iC-jC, ic1, ic2, neqb,eqstring));
                            }
                        }
                    }
                }
                // printf("iC=%d mv.h\n",iC);
                // for (ic1=iC;ic1<nC;ic1++) {
                //     for (ic2=iC;ic2<nC;ic2++) {
                //         printf("h1[ih(%d,%d,%d)=%d]=%1.5f - h2[ih(%d,%d,%d)=%d]=%1.5f = %1.10f\n",
                //                ic1,ic2,iC,ih(nC,ic1,ic2,iC),mv.h1[ih(nC,ic1,ic2,iC)],
                //                ic2,ic1,iC,ih(nC,ic2,ic1,iC),mv.h2[ih(nC,ic2,ic1,iC)],
                //                mv.h1[ih(nC,ic1,ic2,iC)]-mv.h2[ih(nC,ic2,ic1,iC)]);
                //     }
                // }
                // printf("iC=%d hc\n",iC);
                // for (ic1=iC;ic1<nC;ic1++) {
                //     for (ic2=iC;ic2<nC;ic2++) {
                //         printf("hc1[ih(%d,%d,%d)=%d]=%1.5f - hc2[ih(%d,%d,%d)=%d]=%1.5f = %1.10f\n",
                //                ic1,ic2,iC,ih(nC,ic1,ic2,iC),hc1[ih(nC,ic1,ic2,iC)],
                //                ic2,ic1,iC,ih(nC,ic2,ic1,iC),hc2[ih(nC,ic2,ic1,iC)],
                //                hc1[ih(nC,ic1,ic2,iC)]-hc2[ih(nC,ic2,ic1,iC)]);
                //     }
                // }
            }
        }
        /********************************************************/
        // SECTION 3.0: Solve the (c,c,c) corner game
        /********************************************************/
        ic1=iC; ic2=iC; ieqb=0; // Start c=c1=c2=0, only one equilibrium in the corner  endgame
        
        g[iC].pf1[IdOUT]=bne.pf1[IdBNE];
        g[iC].pf2[IdOUT]=bne.pf2[IdBNE];
        g[iC].ic1[IdOUT]=ic1;
        g[iC].ic2[IdOUT]=ic2;
        g[iC].ieqb[IdOUT]=ieqb;
        g[iC].neqb[IdOUT]=1;
        g[iC].c1[IdOUT]=mp.cgrid[ic1];
        g[iC].c2[IdOUT]=mp.cgrid[ic2];
        g[iC].eqbtype[IdOUT]=1;    // 0 = pure strategy equilibrium
        
        // Solve (c,c,c) endgame for firm 1
        g[iC].v10[IdOUT]=(bne.pf1[IdBNE]+mp.df*mv.pc*mv.h1[ih(nC, iC, iC, iC)]+mp.df*(1-mv.pc)*mv.logsumK)/(1-mp.df*(1-mv.pc));
        // Solve (c,c,c) endgame for firm 2
        g[iC].v20[IdOUT]=(bne.pf2[IdBNE]+mp.df*mv.pc*mv.h2[ih(nC, iC, iC, iC)]+mp.df*(1-mv.pc)*mv.logsumK)/(1-mp.df*(1-mv.pc));
        
        g[iC].v11[IdOUT]=g[iC].v10[IdOUT]-mv.kv;
        g[iC].v21[IdOUT]=g[iC].v20[IdOUT]-mv.kv;
        
        g[iC].p1[IdOUT]=f_logit(nF, mp.eta, NULL, g[iC].v11[IdOUT], g[iC].v10[IdOUT]);
        g[iC].p2[IdOUT]=f_logit(nF, mp.eta, NULL, g[iC].v21[IdOUT], g[iC].v20[IdOUT]);
        g[iC].seleqb[id(nC,ic1,ic2,iC,ieqb)]=1;

 		// Expected costs, (c,c,c) corner game
        logsum1=f_lnsum(2, mp.eta, g[iC].v10[IdOUT], g[iC].v11[IdOUT]);		
        logsum2=f_lnsum(2, mp.eta, g[iC].v20[IdOUT], g[iC].v21[IdOUT]);		
		EP11 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v11[IdOUT])+mv.kv
                +mp.df*mv.pc*hc1[ih(nC, iC, iC, iC)];
		EP10 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v10[IdOUT])
                +mp.df*mv.pc*hc1[ih(nC, iC, iC, iC)];
		EP21 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v21[IdOUT])+mv.kv
                +mp.df*mv.pc*hc2[ih(nC, iC, iC, iC)];
		EP20 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v20[IdOUT])
                +mp.df*mv.pc*hc2[ih(nC, iC, iC, iC)];
        p1=g[iC].p1[IdOUT];
        p2=g[iC].p2[IdOUT];

        //firm 1 symbolic solution
        j=ih(nC, ic1, ic2, iC);
        B=mp.df*(1-mv.pc); //used everywhere below!!!!
        AI1=mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v11[IdOUT])+mv.kv
           +mp.df*mv.pc*hc1[j];
        AN1=mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v10[IdOUT])
           +mp.df*mv.pc*hc1[j];

        g[iC].ecv10[IdOUT]=-(AN1 + AI1*B*p1 - AN1*B*p1)/(B - 1);//CN1_ccc
        g[iC].ecv11[IdOUT]=-(AI1 - AI1*B + AN1*B + AI1*B*p1 - AN1*B*p1)/(B - 1);//CI1_ccc

        //firm 1 analytic solution
        // g[iC].ecv11[IdOUT]=(EP11+((mp.df*(1-mv.pc)*(1-p1))/(1-(1-p1)*mp.df*(1-mv.pc)))*EP10)/(1-(mp.df*(1-mv.pc)*p1+((mp.df*(1-mv.pc)*(1-p1)*mp.df*(1-mv.pc)*p1)/(1-(1-p1)*mp.df*(1-mv.pc)))));
        // g[iC].ecv10[IdOUT]=(EP10+mp.df*(1-mv.pc)*p1*g[iC].ecv11[IdOUT])/(1-(1-p1)*mp.df*(1-mv.pc));
        
        //firm 2 symbolic solution
        AI2=mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v21[IdOUT])+mv.kv
           +mp.df*mv.pc*hc2[j];
        AN2=mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v20[IdOUT])
           +mp.df*mv.pc*hc2[j];
// CN2_ccc =
// -(AN2 + AI2*B*p2 - AN2*B*p2)/(B - 1)
// CI2_ccc =
// -(AI2 - AI2*B + AN2*B + AI2*B*p2 - AN2*B*p2)/(B - 1)

        //firm 2 analytic solution
        g[iC].ecv21[IdOUT]=(EP21+((mp.df*(1-mv.pc)*(1-p2))/(1-(1-p2)*mp.df*(1-mv.pc)))*EP20)/(1-(mp.df*(1-mv.pc)*p2+((mp.df*(1-mv.pc)*(1-p2)*mp.df*(1-mv.pc)*p2)/(1-(1-p2)*mp.df*(1-mv.pc)))));
        g[iC].ecv20[IdOUT]=(EP20+mp.df*(1-mv.pc)*p2*g[iC].ecv21[IdOUT])/(1-(1-p2)*mp.df*(1-mv.pc));
        
        /**************************************************************************/
        // SECTION 3.1: solve the (c1,c,c) edge game
        /**************************************************************************/
        ic2=iC; ieqb=0; //only one equilibrium in the edge  endgame
        for (ic1=iC+1;ic1<nC;ic1++) {
            if (((ic1z==ic1) && (ic2z==iC)) || (iC>iCz) || (iCz==0)) {
                // Solve for BN-equaliprium prices and profits
                // Part of formula ??
                // Save output
                g[iC].pf1[IdOUT]=bne.pf1[IdBNE];
                g[iC].pf2[IdOUT]=bne.pf2[IdBNE];
                g[iC].ic1[IdOUT]=ic1;
                g[iC].ic2[IdOUT]=ic2;
                g[iC].ieqb[IdOUT]=ieqb;
                g[iC].neqb[IdOUT]=1;
                g[iC].c1[IdOUT]=mp.cgrid[ic1];
                g[iC].c2[IdOUT]=mp.cgrid[ic2];
                g[iC].eqbtype[IdOUT]=1;
                g[iC].seleqb[IdOUT]=1;
                
                // solve value functions for firm 1
                q=mv.pc*mv.h1[ih(nC, iC, iC, iC)]+(1-mv.pc)*g[iC].v10[id(nC, iC, iC, iC, 0)];
                g[iC].v11[IdOUT]=bne.pf1[IdBNE]-mv.kv+mp.df*q+mp.df*(1-mv.pc)*mv.logsumK;
                g[iC].v10[IdOUT]=f_solvevf(
                        bne.pf1[IdBNE]+mp.df*mv.pc*mv.h1[ih(nC, ic1, iC, iC)],
                        (1-mv.pc)*mp.df,
                        g[iC].v11[id(nC, ic1, iC, iC, 0)],
                        mp.eta, &mp);
                g[iC].p1[IdOUT]=f_logit(nF, mp.eta, NULL, g[iC].v11[IdOUT], g[iC].v10[IdOUT]);
                
                // solve value functions for firm 2
                q=mv.pc*mv.h2[ih(nC, iC, iC, iC)]+(1-mv.pc)*g[iC].v20[id(nC, iC, iC, iC, 0)];
                g[iC].v20[IdOUT]=(bne.pf2[IdBNE]
                        +mp.df*g[iC].p1[IdOUT]*q
                        +mp.df*(1-g[iC].p1[IdOUT])*mv.pc*mv.h2[ih(nC, ic1, iC, iC)]
                        +mp.df*(1-mv.pc)*mv.logsumK
                        )/(1-(1-mv.pc)*mp.df*(1-g[iC].p1[IdOUT]));   // equation (26)
                g[iC].v21[IdOUT]=g[iC].v20[IdOUT]-mv.kv;
                g[iC].p2[IdOUT]=f_logit(nF, mp.eta, NULL, -mv.kv, 0);
                
                 // Expected costs
                logsum1=f_lnsum(2, mp.eta, g[iC].v10[IdOUT], g[iC].v11[IdOUT]);
                logsum2=f_lnsum(2, mp.eta, g[iC].v20[IdOUT], g[iC].v21[IdOUT]);
                EP11 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v11[IdOUT])+mv.kv
                         +mp.df*mv.pc*(hc1[ih(nC, iC, iC, iC)]);
                EP10 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v10[IdOUT])
                         +mp.df*mv.pc*(hc1[ih(nC, ic1, iC, iC)]);
                EP21 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v21[IdOUT])+mv.kv
                         +mp.df*mv.pc*(g[iC].p1[IdOUT]*hc2[ih(nC, iC, iC, iC)]+(1-g[iC].p1[IdOUT])*hc2[ih(nC, ic1, iC, iC)]);
                EP20 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v20[IdOUT])
                         +mp.df*mv.pc*(g[iC].p1[IdOUT]*hc2[ih(nC, iC, iC, iC)]+(1-g[iC].p1[IdOUT])*hc2[ih(nC, ic1, iC, iC)]);
                
                iccc=id(nC, iC, iC, iC, 0);

                //firm1 symbolic solution
                j=ih(nC, ic1, ic2, iC); //index of current state point (without ieqb)
                p1=g[iC].p1[IdOUT];
                p2=g[iC].p2[IdOUT];
                AN1=mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v10[IdOUT])
                   +mp.df*mv.pc*hc1[j];

                g[iC].ecv11[IdOUT] =mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v11[IdOUT])+mv.kv
                                   +mp.df*mv.pc*(hc1[ih(nC, iC, iC, iC)])
                                   +B*(g[iC].p1[iccc]*g[iC].ecv11[iccc]  +(1-g[iC].p1[iccc])*g[iC].ecv10[iccc]); //CI1_c1cc
                g[iC].ecv10[IdOUT] =(AN1 + g[iC].ecv11[IdOUT]*B*p1)/(B*p1 - B + 1); //CN1_c1cc
 
                //firm1 analytic solution
                // g[iC].ecv11[IdOUT] =EP11    +mp.df*(1-mv.pc)*(  g[iC].p1[iccc])*g[iC].ecv11[iccc]
                //                             +mp.df*(1-mv.pc)*(1-g[iC].p1[iccc])*g[iC].ecv10[iccc];

                // g[iC].ecv10[IdOUT] =(EP10  + mp.df*(1-mv.pc)*g[iC].p1[IdOUT]*g[iC].ecv11[IdOUT])
                //                             /(1-mp.df*(1-mv.pc)*(1-g[iC].p1[IdOUT]));


                //firm 2 analytic solution
                f1=(mp.df*(1-mv.pc)*(1-g[iC].p1[IdOUT])*g[iC].p2[IdOUT])/(1-mp.df*(1-mv.pc)*(1-g[iC].p1[IdOUT])*g[iC].p2[IdOUT]);
                g[iC].ecv20[IdOUT]=(EP20+f1*EP20+(mp.df*(1-mv.pc)*g[iC].p1[IdOUT]+f1)*g[iC].p2[iccc]*g[iC].ecv21[iccc]+f1*(mp.df*(1-mv.pc)*g[iC].p1[IdOUT]*(1-g[iC].p2[iccc])*g[iC].ecv20[iccc]))
                                /(1-(1+f1)*mp.df*(1-mv.pc)*(1-g[iC].p1[IdOUT])*(1-g[iC].p2[IdOUT]));
                
                g[iC].ecv21[IdOUT]=(EP21+g[iC].p2[IdOUT]*g[iC].ecv21[iccc]+mp.df*(1-mv.pc)*g[iC].p1[IdOUT]*(1-g[iC].p2[IdOUT])*g[iC].ecv20[iccc]+mp.df*(1-mv.pc)*(1-g[iC].p1[IdOUT])*(1-g[iC].p2[IdOUT])*g[iC].ecv20[IdOUT])
                                /(1-mp.df*(1-mv.pc)*(1-g[iC].p1[IdOUT])*g[iC].p2[IdOUT]);
                
            }
        }
        
        /**************************************************************************/
        // SECTION 3.2: solve the (c,c2,c) edge game
        /**************************************************************************/
        ic1=iC; ieqb=0; //only one equilibrium in the edge endgame
        for (ic2=iC+1;ic2<nC;ic2++) {
            if (((ic2z==ic2) && (ic1z==iC)) || (iC>iCz) || (iCz==0)) {
                
                // Save output
                g[iC].pf1[IdOUT]=bne.pf1[IdBNE];
                g[iC].pf2[IdOUT]=bne.pf2[IdBNE];
                g[iC].ic1[IdOUT]=ic1;
                g[iC].ic2[IdOUT]=ic2;
                g[iC].ieqb[IdOUT]=ieqb;
                g[iC].neqb[IdOUT]=1;
                g[iC].c1[IdOUT]=mp.cgrid[ic1];
                g[iC].c2[IdOUT]=mp.cgrid[ic2];
                g[iC].eqbtype[IdOUT]=1;
                g[iC].seleqb[IdOUT]=1;
                
                // solve value functions for firm 2
                q=mv.pc*mv.h2[ih(nC, iC, iC, iC)]+(1-mv.pc)*g[iC].v20[id(nC, iC, iC, iC, 0)];
                g[iC].v21[IdOUT]=bne.pf2[IdBNE]-mv.kv+mp.df*q+mp.df*(1-mv.pc)*mv.logsumK;
                g[iC].v20[IdOUT]=f_solvevf(
                        bne.pf2[IdBNE]+mp.df*mv.pc*mv.h2[ih(nC, iC, ic2, iC)],
                        (1-mv.pc)*mp.df,
                        g[iC].v21[IdOUT],
                        mp.eta, &mp);
                g[iC].p2[IdOUT]=f_logit(nF, mp.eta, NULL, g[iC].v21[IdOUT], g[iC].v20[IdOUT]);
                
                // solve value functions for firm 1
                q=mv.pc*mv.h1[ih(nC, iC, iC, iC)]+(1-mv.pc)*g[iC].v10[id(nC, iC, iC, iC, 0)];
                g[iC].v10[IdOUT]=(
                        bne.pf1[IdBNE]
                        +mp.df*g[iC].p2[id(nC,iC,ic2,iC,0)]*q
                        +mp.df*(1-g[iC].p2[id(nC,iC,ic2,iC,0)])*mv.pc*mv.h1[ih(nC, iC, ic2, iC)]
                        +mp.df*(1-mv.pc)*mv.logsumK
                        )/(1-(1-mv.pc)*mp.df*(1-g[iC].p2[id(nC,iC,ic2,iC,0)]));   // equation (26)
                g[iC].v11[IdOUT]=g[iC].v10[IdOUT]-mv.kv;
                g[iC].p1[IdOUT]=f_logit(nF, mp.eta, NULL, g[iC].v11[IdOUT], g[iC].v10[IdOUT]);
                
                // Expected costs
                logsum1=f_lnsum(2, mp.eta, g[iC].v10[IdOUT], g[iC].v11[IdOUT]);
                logsum2=f_lnsum(2, mp.eta, g[iC].v20[IdOUT], g[iC].v21[IdOUT]);
                EP11 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v11[IdOUT])+mv.kv
                                         +mp.df*mv.pc*(g[iC].p2[IdOUT]*hc1[ih(nC, iC, iC, iC)]+(1-g[iC].p2[IdOUT])*hc1[ih(nC, iC, ic2, iC)]);
                EP10 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v10[IdOUT])
                                         +mp.df*mv.pc*(g[iC].p2[IdOUT]*hc1[ih(nC, iC, iC, iC)]+(1-g[iC].p2[IdOUT])*hc1[ih(nC, iC, ic2, iC)]);
                EP21 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v21[IdOUT])+mv.kv
                                         +mp.df*mv.pc*(hc2[ih(nC, iC, iC, iC)]);
                EP20 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v20[IdOUT])
                                         +mp.df*mv.pc*(hc2[ih(nC, iC, ic2, iC)]);
                
                iccc=id(nC, iC, iC, iC, 0);

                //firm1 symbolic solution
                j=ih(nC, ic1, ic2, iC); //index of current state point (without ieqb)
                p1=g[iC].p1[IdOUT];
                p2=g[iC].p2[IdOUT];
                AN1=mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v10[IdOUT])
                   +mp.df*mv.pc*(g[iC].p2[IdOUT]*hc1[ih(nC,iC,iC,iC)]+(1-g[iC].p2[IdOUT])*hc1[j])
                   +B*p2*(g[iC].p1[iccc]*g[iC].ecv11[iccc]+(1-g[iC].p1[iccc])*g[iC].ecv10[iccc]);
                AI1=mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v11[IdOUT])+mv.kv
                   +mp.df*mv.pc*(g[iC].p2[IdOUT]*hc1[ih(nC,iC,iC,iC)]+(1-g[iC].p2[IdOUT])*hc1[j])
                  + B*p2*(g[iC].p1[iccc]*g[iC].ecv11[iccc]+(1-g[iC].p1[iccc])*g[iC].ecv10[iccc]);
                 
                g[iC].ecv10[IdOUT]=(AN1 + AI1*B*p1 - AN1*B*p1 - AI1*B*p1*p2 + AN1*B*p1*p2)/(B*p2 - B + 1);//CN1_cc2c
                g[iC].ecv11[IdOUT]=(AI1 - AI1*B + AN1*B + AI1*B*p1 + AI1*B*p2 - AN1*B*p1 - AN1*B*p2 - AI1*B*p1*p2 + AN1*B*p1*p2)/(B*p2 - B + 1);//CI1_cc2c

                //firm1 analytic solution
                // f1=(mp.df*(1-mv.pc)*(1-g[iC].p2[IdOUT])*g[iC].p1[IdOUT])/(1-mp.df*(1-mv.pc)*(1-g[iC].p2[IdOUT])*g[iC].p1[IdOUT]);

                // g[iC].ecv10[IdOUT]=(EP10+f1*EP10+(mp.df*(1-mv.pc)*g[iC].p2[IdOUT]+f1)*g[iC].p1[iccc]*g[iC].ecv11[iccc]+f1*(mp.df*(1-mv.pc)*g[iC].p2[IdOUT]*(1-g[iC].p1[iccc])*g[iC].ecv10[iccc]))
                //                     /(1-(1+f1)*mp.df*(1-mv.pc)*(1-g[iC].p2[IdOUT])*(1-g[iC].p1[IdOUT]));
                
                // g[iC].ecv11[IdOUT]=(EP11+g[iC].p1[IdOUT]*g[iC].ecv11[iccc]+mp.df*(1-mv.pc)*g[iC].p2[IdOUT]*(1-g[iC].p1[IdOUT])*g[iC].ecv10[iccc]+mp.df*(1-mv.pc)*(1-g[iC].p2[IdOUT])*(1-g[iC].p1[IdOUT])*g[iC].ecv10[IdOUT])
                //                     /(1-mp.df*(1-mv.pc)*(1-g[iC].p2[IdOUT])*g[iC].p1[IdOUT]);

                //firm 2 analytic solution
                g[iC].ecv21[IdOUT] =EP21    +mp.df*(1-mv.pc)*(  g[iC].p2[iccc])*g[iC].ecv21[iccc]
                                            +mp.df*(1-mv.pc)*(1-g[iC].p2[iccc])*g[iC].ecv20[iccc];

                g[iC].ecv20[IdOUT] =(EP20  + mp.df*(1-mv.pc)*g[iC].p2[IdOUT]*g[iC].ecv21[IdOUT])
                                            /(1-mp.df*(1-mv.pc)*(1-g[iC].p2[IdOUT]));
           }
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
                if (((ic1z==ic1) && (ic2z==ic2)) || (iC>iCz) || (iCz==0)) {
                    j=(ic1-iC)*(nC-iC)+(ic2-iC);    // index used for storage in h
                    i1=id(nC, ic1, iC, iC, 0); //c1 c c
                    i2=id(nC, iC, ic2, iC, 0); //c c2 c
                    
                    // populate brc sturcture
                    brc.pc=mv.pc;        // XXX appears in two structures
                    brc.kv=mv.kv;        // XXX appears in two structures
                    brc.pf[0]=bne.pf1[IdBNE];
                    brc.pf[1]=bne.pf2[IdBNE];
                    
                    brc.x0[0]= mp.df*(brc.pc*mv.h1[ih(nC,iC,iC,iC)]     +    (1-brc.pc)*g[iC].v10[id(nC,iC,iC,iC,0)]);
                    brc.x1[0]= mp.df*(brc.pc*mv.h1[ih(nC,iC,ic2,iC)]    +    (1-brc.pc)*g[iC].v10[id(nC,iC,ic2,iC,0)]);
                    
                    brc.x0[1]= mp.df*(brc.pc*mv.h2[ih(nC,iC,iC,iC)]     +    (1-brc.pc)*g[iC].v20[id(nC,iC,iC,iC,0)]);
                    brc.x1[1]= mp.df*(brc.pc*mv.h2[ih(nC,ic1,iC,iC)]    +    (1-brc.pc)*g[iC].v20[id(nC,ic1,iC,iC,0)]);
                    
                    brc.A0[0]=        bne.pf1[IdBNE]+mp.df*mv.pc*mv.h1[ih(nC, iC, iC, iC)]+(1-mv.pc)*f_lnsum(nF, mp.eta, g[iC].v10[i1], g[iC].v11[i1]);
                    brc.A1[0]=        mp.df*mv.pc*mv.h1[ih(nC,ic1,ic2,iC)];
                    
                    brc.A0[1]=        bne.pf2[IdBNE]+mp.df*mv.pc*mv.h2[ih(nC, iC, iC, iC)]+(1-mv.pc)*f_lnsum(nF, mp.eta, g[iC].v20[i2], g[iC].v21[i2]);
                    brc.A1[1]=        mp.df*mv.pc*mv.h2[ih(nC,ic1,ic2,iC)];
                    
                    brc.B[1]=mp.df*(1-brc.pc);
                    brc.B[0]=mp.df*(1-brc.pc);
                    
                    //initialize the number of equilibria to be stored in the g structure
                    g[iC].neqb[id(nC,ic1,ic2,iC,0)]=0;
                                        
                    //find equilibria
                    if ((mp.eta>0) && (mp.analytical==false)){
                         s_ffxp((f_br2p1_sm), &mp, &brc, &g[iC], ic1, ic2, iC, nC);
//                         s_ffxp((f_br1p2_sm), &mp, &brc, &g[iC], ic1, ic2, iC, nC);
                    }
                    else {
                        analytical_solutions (ic1,ic2,iC,mp,mv,bne,g,br,eqstring);
                    }

                    
                    //number of equilibria found
                    //neqb=f_neqb(nC, ic1, ic2, iC, g); XXXX
                    neqb=(int)g[iC].neqb[id(nC,ic1,ic2,iC,0)];
#ifdef PRINTnrstagempe
                    if (neqb==PRINTnrstagempe) {
                        neqstr=nC*(nC+1)*(2*nC+1)/6;
                        printf("%d stage MPE: iC=%d, ic1=%d, ic2=%d ESS= ",PRINTnrstagempe,iC+1,ic1+1,ic2+1);
                        printeqb(eqstring,neqstr,lexindex1(eqstring,neqstr),NULL);
                    };
#endif
                    
                    // Save indices and costs (remaining elements of g are stored inside f_br2p1_am
                    for (ieqb=0;ieqb<neqb;ieqb++) {
                        g[iC].ic1[IdOUT]=ic1;               //0
                        g[iC].ic2[IdOUT]=ic2;               //1
                        g[iC].ieqb[IdOUT]=ieqb;             //2
                        g[iC].c1[IdOUT]=mp.cgrid[ic1];      //3
                        g[iC].c2[IdOUT]=mp.cgrid[ic2];      //4
                        // 5-11: eqbtype, v10,v11,b20, v21, p1, p2: stored inside f_br2p1_am
                        g[iC].seleqb[id(nC,ic1,ic2,iC,ieqb)]=(ieqb==f_ers(&mp, iC, ic1, ic2, neqb,eqstring)?1.0:0.0);  //12
                        g[iC].pf1[IdOUT]=bne.pf1[IdBNE];           //13
                        g[iC].pf2[IdOUT]=bne.pf2[IdBNE];           //14
                        
                        
                        //calculate analytical solutions and save to proper output structures
                        // Expected costs
                        logsum1=f_lnsum(2, mp.eta, g[iC].v10[IdOUT], g[iC].v11[IdOUT]);
                        logsum2=f_lnsum(2, mp.eta, g[iC].v20[IdOUT], g[iC].v21[IdOUT]);
                        EP11 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v11[IdOUT])+mv.kv
                                                 +mp.df*mv.pc*(g[iC].p2[IdOUT]*hc1[ih(nC, iC , iC, iC)]+(1-g[iC].p2[IdOUT])*hc1[ih(nC, iC , ic2, iC)]); ///BUG HERE
                                              // +mp.df*mv.pc*(g[iC].p2[IdOUT]*hc1[ih(nC, ic1, iC, iC)]+(1-g[iC].p2[IdOUT])*hc1[ih(nC, ic1, ic2, iC)]); ///BUG HERE
                        EP21 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v21[IdOUT])+mv.kv
                                                +mp.df*mv.pc*(g[iC].p1[IdOUT]*hc2[ih(nC, iC, iC, iC)]+(1-g[iC].p1[IdOUT])*hc2[ih(nC, ic1, iC, iC)]);
                                                // +mp.df*mv.pc*(g[iC].p1[IdOUT]*hc2[ih(nC, iC, ic2, iC)]+(1-g[iC].p1[IdOUT])*hc2[ih(nC, ic1, ic2, iC)]); ///BUG HERE
                        EP10 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v10[IdOUT])
                                                +mp.df*mv.pc*(g[iC].p2[IdOUT]*hc1[ih(nC, ic1, iC, iC)]+(1-g[iC].p2[IdOUT])*hc1[ih(nC, ic1, ic2, iC)]);
                        EP20 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v20[IdOUT])
                                                +mp.df*mv.pc*(g[iC].p1[IdOUT]*hc2[ih(nC, iC, ic2, iC)]+(1-g[iC].p1[IdOUT])*hc2[ih(nC, ic1, ic2, iC)]);

                        iccc=id(nC, iC, iC, iC, 0); //0 for first and only equilibrium at corners and edges
                        i1=id(nC, ic1, iC, iC, 0);
                        i2=id(nC, iC, ic2, iC, 0);
                        
                        //firm1 symbolic solution
                        j=ih(nC, ic1, ic2, iC); //index of current state point (without ieqb)
                        p1=g[iC].p1[IdOUT];
                        p2=g[iC].p2[IdOUT];

                        g[iC].ecv11[IdOUT]=mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v11[IdOUT])+mv.kv
                                          +mp.df*mv.pc*(p2 *hc1[ih(nC, iC, iC, iC)]
                                                    +(1-p2)*hc1[ih(nC, iC, ic2, iC)])
                                          +B*p2*    (g[iC].p1[iccc]*g[iC].ecv11[iccc]+(1-g[iC].p1[iccc])*g[iC].ecv10[iccc])
                                          +B*(1-p2)*(g[iC].p1[i2]  *g[iC].ecv11[i2]  +(1-g[iC].p1[i2])  *g[iC].ecv10[i2]); //CI1_c1c2c

                        AN1=               mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v10[IdOUT])
                                          +mp.df*mv.pc*(p2*hc1[ih(nC, ic1, iC, iC)]+(1-p2)*hc1[j])
                                          +B*p2*(g[iC].p1[i1]*g[iC].ecv11[i1]+(1-g[iC].p1[i1])*g[iC].ecv10[i1]);

                        g[iC].ecv10[IdOUT]=(AN1 + g[iC].ecv11[IdOUT]*B*p1 - g[iC].ecv11[IdOUT]*B*p1*p2)/(B*p1 - B + B*p2 - B*p1*p2 + 1);//CN1_c1c2c
                         
                         
                        //firm1 analytic solution
                        // g[iC].ecv11[IdOUT]=EP11 + mp.df*(1-mv.pc)*g[iC].p2[IdOUT]     *(g[iC].p1[iccc]*g[iC].ecv11[iccc]+(1-g[iC].p1[iccc])*g[iC].ecv10[iccc])
                        // + mp.df*(1-mv.pc)*(1-g[iC].p2[IdOUT]) *(g[iC].p1[i2]*g[iC].ecv11[i2]   +(1-g[iC].p1[i2])  *g[iC].ecv10[i2]);

                        // g[iC].ecv10[IdOUT]=(EP10+mp.df*(1-mv.pc)*   g[iC].p2[IdOUT] *(g[iC].p1[i1]* g[iC].ecv11[i1]
                        //                                                               +(1-g[iC].p1[i1])*g[iC].ecv10[i1])
                        //                         +mp.df*(1-mv.pc)*(1-g[iC].p2[IdOUT]) *g[iC].p1[IdOUT]*g[iC].ecv11[IdOUT])
                        //                     /(1-mp.df*(1-mv.pc)*(1-g[iC].p2[IdOUT])*(1-g[iC].p1[IdOUT]));
                                                
                        //firm2 symbolic solution
                        // g[iC].ecv21[IdOUT]=mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v21[IdOUT])+mv.kv
                        //                   +mp.df*mv.pc*(p1*hc2[ih(nC, iC, iC, iC)]+(1-p1)*hc2[ih(nC, ic1, iC, iC)])
                        //                   +B*p1    *(g[iC].p2[iccc]*g[iC].ecv21[iccc]+(1-g[iC].p2[iccc])*g[iC].ecv20[iccc])
                        //                   +B*(1-p1)*(g[iC].p2[i1]*g[iC].ecv21[i1]   +(1-g[iC].p2[i1])  *g[iC].ecv20[i1]);//CI2_c1c2c
                        
                        // AN2=               mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v20[IdOUT])
                        //                   +mp.df*mv.pc*(p1*hc2[ih(nC, iC, ic2, iC)]+(1-p1)*hc2[j])
                        //                   +B*p1*(g[iC].p2[i2]*g[iC].ecv21[i2]+(1-g[iC].p2[i2])*g[iC].ecv20[i2]);
                         
                        // g[iC].ecv20[IdOUT]=(AN2 + g[iC].ecv21[IdOUT]*B*p2 - g[iC].ecv21[IdOUT]*B*p1*p2)/(B*p1 - B + B*p2 - B*p1*p2 + 1);//CN2_c1c2c

                        //firm2 analytic solution
                        g[iC].ecv21[IdOUT]=EP21 + mp.df*(1-mv.pc)*g[iC].p1[IdOUT]     *(g[iC].p2[iccc]*g[iC].ecv21[iccc]+(1-g[iC].p2[iccc])*g[iC].ecv20[iccc])
                        + mp.df*(1-mv.pc)*(1-g[iC].p1[IdOUT]) *(g[iC].p2[i1]*g[iC].ecv21[i1]   +(1-g[iC].p2[i1])  *g[iC].ecv20[i1]);
                        
                        g[iC].ecv20[IdOUT]=(EP20+mp.df*(1-mv.pc)*g[iC].p1[IdOUT] *(g[iC].p2[i2]*g[iC].ecv21[i2]+(1-g[iC].p2[i2])*g[iC].ecv20[i2])
                                             +mp.df*(1-mv.pc)*(1-g[iC].p1[IdOUT])*g[iC].p2[IdOUT]*g[iC].ecv21[IdOUT] )
                                            /(1-mp.df*(1-mv.pc)*(1-g[iC].p1[IdOUT])*(1-g[iC].p2[IdOUT]));

                    }
                    

                            
                            
//                             g[iC].ecv10[IdOUT]=(EP10+f1*EP10+(mp.df*g[iC].p2[IdOUT]+f1)*g[iC].p1[iccc]*g[iC].ecv11[iccc]+f1*(mp.df*g[iC].p2[IdOUT]*(1-g[iC].p1[iccc])*g[iC].ecv10[iccc]))
//                     /(1-(1+f1)*mp.df*(1-g[iC].p2[IdOUT])*(1-g[iC].p1[IdOUT]));
                    
                    
                    


                    
                }
            }
        }
//     }  // end parallel execution

        //compute totals for the whole layer
        ieqb=0; //only one equilibrium in the edge endgame
        for (ic1=iC;ic1<nC;ic1++) {
            for (ic2=iC;ic2<nC;ic2++) {
                if (((ic1z==ic1) && (ic2z==ic2)) || (iC>iCz) || (iCz==0)) {
                    if (ic1>iC && ic2>iC) neqb=(int)g[iC].neqb[id(nC,ic1,ic2,iC,0)]; //interior
                    else neqb=1; //corner and edges 
                    for (ieqb=0;ieqb<neqb;ieqb++) {
                        //social surplus = total revenue - total cost (conditional on who has the right of move)
                        g[iC].ecm1[IdOUT]=mp.dt*MAX(mp.cgrid[ic1],mp.cgrid[ic2])/(1-mp.df) - (g[iC].p1[IdOUT]*g[iC].ecv11[IdOUT] + (1-g[iC].p1[IdOUT])*g[iC].ecv10[IdOUT]
                                                                                             +g[iC].p2[IdOUT]*g[iC].ecv21[IdOUT] + (1-g[iC].p2[IdOUT])*g[iC].ecv20[IdOUT]);
                        //printf("neqb=%d, ieqb=%d,  p1=%1.4f,  ecv11=%1.4f,  ecv10=%1.4f,   ecm1=%1.4f,  ecm2=%1.4f\n",neqb ,ieqb, g[iC].p1[IdOUT], g[iC].ecv11[IdOUT], g[iC].ecv10[IdOUT],g[iC].ecm1[IdOUT],g[iC].ecm2[IdOUT] );
                    }
                }
            }
        }
    } //continue loop over iC

    
    free(mv.h1);
    free(mv.h2);
    free(hc1);
    free(hc2);
    free(mv.eqbinfo);
//printf("Run %d\n",COUNT++);
}
