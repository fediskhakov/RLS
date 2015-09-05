// Method for looking for equilibria: EQBMETHOD_GRID, EQBMETHOD_VF
#define EQBMETHOD_GRID

// Switch for drawing equilibrium graph EQB_VISUAL
// #define EQB_VISUAL

// #define OTHER_INPUTS ,Gamestruct *g,MPstruct *mp,MVstruct *mv,int ic1,int ic2,int iC,int nC
// #define OTHER_INPUTS_CLEAN ,g,mp,mv,ic1,ic2,iC,nC

// void s_FindEqb (double (*fun)(double p1, int eqnum, int eqtype OTHER_INPUTS) OTHER_INPUTS) {

// double f_brp(double p1, int ieqb, int eqbtype, Gamestruct *g, MPstruct *mp, MVstruct *mv, int ic1, int ic2, int iC, int nC) {
// f_br2p1_am(double p1, int id, MPstruct* mp, double* pf1, double* pf2, Gamestruct* gc) {


double f_br2p1_am(double p1, int ieqb, int eqbtype, MPstruct* mp, Brc_amstruct* brc, Gamestruct* gc, int ic1, int ic2, int iC, int nC);
double f_brp1_am(double p2, MPstruct* mp, Brc_amstruct* brc);
double f_brp2_am(double p1, MPstruct* mp, Brc_amstruct* brc);
void s_plot_bpr(MPstruct* mp, Brc_amstruct* brc, int ic1, int ic2, int iC) ;

void leap_am(MPstruct mp, Gamestruct* g, Bnestruct bne, Brstruct* br, int* eqstring, int ic1z, int ic2z, int iCz) {
    /********************************************************/
// SECTION 1: DECLARATIONS
    /********************************************************/
    // Params and Modparams
    int nC; //Main parameter
    MVstruct mv;
    Brc_amstruct brc;
    
    //double v, v1, r, p1;
    int i, j, iC, jC, ic1, ic2, ieqb, neqb, i1, i2;
    //double q; //temps for c>0
    //double v, v1, r, p1;
    //int i, j, iF, iI, iC, ic1, ic2, ieqb, it, ip, neqb;
    double f1, f2, A,B,C, EP10, EP11, EP20, EP21,EPx10, EPx1, EPx2, logsum1, logsum2, p1, p2;
    // XXX introduce as parameter
    double f11, f22, fp1, fp2;
    double A1, A2, B1, B2, C1, C2;   // auxiliary parameters used for readability in cost recursion in c1 c2 c game
    double P1, P2, X, CI1, CI2;
    double B11, B12, B21, B22, AN1v, AN1x, AN2v, AN2x;
    double *h11; double *h12; double *h21; double *h22;
    double *hc11; double *hc12; double *hc21; double *hc22;
    int jp;

    mp.tpm11=mp.tpm[0]; mp.tpm21=mp.tpm[2];   
    mp.tpm12=mp.tpm[1]; mp.tpm22=mp.tpm[3];

    
    nC        =  mp.nC;               // very much used => for easy referencing
    
    /********************************************************/
// SECTION 2: MEMORY ALLOCATION
    /********************************************************/
    h11      = (double *) calloc(nC*nC, sizeof(double));     // temp storage of formula parts, precomputing, see c>0
    h12      = (double *) calloc(nC*nC, sizeof(double));     // temp storage of formula parts, precomputing, see c>0
    h21      = (double *) calloc(nC*nC, sizeof(double));     // temp storage of formula parts, precomputing, see c>0
    h22      = (double *) calloc(nC*nC, sizeof(double));     // temp storage of formula parts, precomputing, see c>0
    hc11      = (double *) calloc(nC*nC, sizeof(double));     // temp storage of formula parts, precomputing, see c>0
    hc12      = (double *) calloc(nC*nC, sizeof(double));     // temp storage of formula parts, precomputing, see c>0
    hc21      = (double *) calloc(nC*nC, sizeof(double));     // temp storage of formula parts, precomputing, see c>0
    hc22      = (double *) calloc(nC*nC, sizeof(double));     // temp storage of formula parts, precomputing, see c>0
    mv.eqbinfo = (double *) calloc(5*MAXEQBINFO, sizeof(double));     // storage for found equilibria in findeqb
    
// printf("Solving alternating move game\n");
    
    
//    1) If we recompute an interior point, do not recompute edges and corner                   --- if (ic1z!=iCz) && (ic2z!=iCz)
//    2) If we recompute a point on (c1,c,c) edge, do not recompute corner and any other poinst ad eges --- if (ic2z!=iCz) 
//       If we recompute a point on (c,c2,c) edge, do not recompute corner and any other poinst ad eges ---  (ic1z!=iCz)
//    2) If we recompute an edge, do not recompute corner and any other poinst ad eges          --- if (ic1z!=iCz) || (ic2z!=iCz)
//    3) otherwise both edges and corner                                                        --- if (ic1z==iCz) && (ic2z==iCz)  
    
   
        
    /********************************************************/
// SECTION 3: LOOP OVER iC
    /********************************************************/
    for (iC=iCz;iC<nC;iC++) { 
        
        /********************************************************/
        // Initialize for iC
        /********************************************************/
        g[iC].c[0]=mp.cgrid[iC]; // save value of c in output
        
        mv.pc=f_pti(iC,&mp);     // probability of a technological improvement. Initially improvements are assumed to be
        // Precompuation of variables used in cost recursions, that does not depend on state variables 
        f11=(mp.df*(1-mv.pc))/(1-mp.df*(1-mv.pc)*mp.tpm11);
        f22=(mp.df*(1-mv.pc))/(1-mp.df*(1-mv.pc)*mp.tpm22);
        

        // 'incremental', that is the jump, if it occurs, is to c'=cgrid(i-1) with probability pti
        // for iC=0 pti=0!
        mv.kv=f_kf(iC, &mp);     // Investment costs
        mv.logsumK=f_lnsum(2, mp.eta, 0, -mv.kv);        // Update mv.kv using investments costs, K(c)
        
        // precompute h(c1,c2,c) first element of H(c1,c2,c);
//        omp_set_num_threads(NUM_THREADS);
//        #pragma omp parallel shared(mv, g, mp,iC,  nC) private(ic1, ic2, neqb, i, j)
//               printf("There are %d threads\n",omp_get_num_threads());
//        #pragma omp for
        for (ic1=iC;ic1<nC;ic1++) {
            for (ic2=iC;ic2<nC;ic2++) {
                j=(ic1-iC)*(nC-iC)+(ic2-iC);    // index used for storage in h
                h11[j]=0;
                h12[j]=0;
                h21[j]=0;
                h22[j]=0;

                hc11[j]=0;
                hc12[j]=0;
                hc21[j]=0;
                hc22[j]=0;
                for (jC=1;jC<=iC;jC++) {
                    jp=(iC-jC)*nC+iC;
                    if (mp.pti[jp]>0) {
                        neqb=(int)g[iC-jC].neqb[id(nC,ic1,ic2,iC-jC,0)];
                        i=id(nC, ic1, ic2, iC-jC, f_ers(&mp, iC-jC, ic1, ic2, neqb,eqstring));    // index for selected equilibrium
                        
                        h11[j]=h11[j]+mv.pc*mp.df*mp.pti[jp]*(mp.tpm11*f_lnsum(nF, mp.eta, g[iC-jC].v11[i], g[iC-jC].v10[i])
                        +mp.tpm21*(g[iC-jC].p2[i]*g[iC-jC].x11[i] + (1-g[iC-jC].p2[i])*g[iC-jC].x10[i]));
                        
                        h12[j]=h12[j]+mv.pc*mp.df*mp.pti[jp]*(mp.tpm12*f_lnsum(nF, mp.eta, g[iC-jC].v11[i], g[iC-jC].v10[i])
                        +mp.tpm22*(g[iC-jC].p2[i]*g[iC-jC].x11[i] + (1-g[iC-jC].p2[i])*g[iC-jC].x10[i]));
                        
                        h21[j]=h21[j]+mv.pc*mp.df*mp.pti[jp]*(mp.tpm21*f_lnsum(nF, mp.eta, g[iC-jC].v21[i], g[iC-jC].v20[i])
                        +mp.tpm11*(g[iC-jC].p1[i]*g[iC-jC].x21[i] + (1-g[iC-jC].p1[i])*g[iC-jC].x20[i]));
                        
                        h22[j]=h22[j]+mv.pc*mp.df*mp.pti[jp]*(mp.tpm22*f_lnsum(nF, mp.eta, g[iC-jC].v21[i], g[iC-jC].v20[i])
                        +mp.tpm12*(g[iC-jC].p1[i]*g[iC-jC].x21[i] + (1-g[iC-jC].p1[i])*g[iC-jC].x20[i]));
                        
                        // hcjm is the conditional expectation of the discounted costs for firm j, if firm m has the right to move
                        // hc1m: Discounted costs of firm 1, given that firm m=1,2 moves and technology improves
                        hc11[j]=hc11[j]+mp.df*mv.pc*mp.pti[jp]*(g[iC-jC].p1[i]*g[iC-jC].ecv11[i]+(1-g[iC-jC].p1[i])*g[iC-jC].ecv10[i]);  // firm 1 moves
                        hc12[j]=hc12[j]+mp.df*mv.pc*mp.pti[jp]*(g[iC-jC].p2[i]*g[iC-jC].ecx11[i]+(1-g[iC-jC].p2[i])*g[iC-jC].ecx10[i]);  // firm 2 moves
                        
                        // hc2m: Discounted costs of firm 2, given that firm m=1,2 moves and technology improves
                        hc21[j]=hc21[j]+mp.df*mv.pc*mp.pti[jp]*(g[iC-jC].p1[i]*g[iC-jC].ecx21[i]+(1-g[iC-jC].p1[i])*g[iC-jC].ecx20[i]);  // firm 1 moves
                        hc22[j]=hc22[j]+mp.df*mv.pc*mp.pti[jp]*(g[iC-jC].p2[i]*g[iC-jC].ecv21[i]+(1-g[iC-jC].p2[i])*g[iC-jC].ecv20[i]);  // firm 2 moves
                    }
                }
            }
        }
        
//     printf("kv=%f\n",mv.kv);
        

        
        /********************************************************/
        // SECTION 3.0: Solve the (c,c,c) corner game
        /********************************************************/
        ic1=iC; ic2=iC; ieqb=0; // Start c=c1=c2=0, only one equilibrium in the corner  endgame
        j=(ic1-iC)*(nC-iC)+(ic2-iC);    // index used for storage in h
        
        g[iC].ic1[IdOUT]=ic1;
        g[iC].ic2[IdOUT]=ic2;
        g[iC].c1[IdOUT]=mp.cgrid[ic1];
        g[iC].c2[IdOUT]=mp.cgrid[ic2];
        g[iC].pf1[IdOUT]=bne.pf1[IdBNE];
        g[iC].pf2[IdOUT]=bne.pf2[IdBNE];
        g[iC].ieqb[IdOUT]=ieqb;
        g[iC].neqb[IdOUT]=1;
        g[iC].eqbtype[IdOUT]=1;   // 0 = pure strategy equilibrium
        g[iC].seleqb[IdOUT]=1;
        
        // Precompute coefficients
        B11=  (mp.df*(1-mv.pc))*mp.tpm11;
        B21=  (mp.df*(1-mv.pc))*mp.tpm21;
        B12=  (mp.df*(1-mv.pc))*mp.tpm12;
        B22=  (mp.df*(1-mv.pc))*mp.tpm22;
        AN1v = bne.pf1[IdBNE] + h11[j]+B11*mv.logsumK;
        AN1x=  bne.pf1[IdBNE] + h12[j]+B12*mv.logsumK;
        AN2v = bne.pf2[IdBNE] + h22[j]+B22*mv.logsumK;
        AN2x=  bne.pf2[IdBNE] + h21[j]+B21*mv.logsumK;

        // value functions for firm 1
        g[iC].v10[IdOUT]=-(AN1v - AN1v*B22 + AN1x*B21)/(B11 + B22 - B11*B22 + B12*B21 - 1); 
        g[iC].x10[IdOUT]=-(AN1x + AN1v*B12 - AN1x*B11)/(B11 + B22 - B11*B22 + B12*B21 - 1);
        g[iC].v11[IdOUT]=g[iC].v10[IdOUT]-mv.kv;
        g[iC].x11[IdOUT]=g[iC].x10[IdOUT];

        // value functions for firm 2
        g[iC].v20[IdOUT]=-(AN2v - AN2v*B22 + AN2x*B21)/(B22 + B11 - B22*B11 + B21*B12 - 1); 
        g[iC].x20[IdOUT]=-(AN2x + AN2v*B12 - AN2x*B22)/(B22 + B11 - B22*B11 + B21*B12 - 1);
        g[iC].v21[IdOUT]=g[iC].v20[IdOUT]-mv.kv;
        g[iC].x21[IdOUT]=g[iC].x20[IdOUT];


        // value functions for firm 1
        // f1=((1-mv.pc)*mp.df*mp.tpm21)/(1-mp.df*(1-mv.pc)*(1-mp.tpm12));
        // g[iC].v10[IdOUT]=((bne.pf1[IdBNE]+h11[j])*(1+f1)+mp.df*(1-mv.pc)*(mp.tpm11+f1*mp.tpm12)*mv.logsumK)/(1-mp.df*(1-mv.pc)*(mp.tpm11+f1*mp.tpm12));
        // g[iC].v11[IdOUT]=g[iC].v10[IdOUT]-mv.kv;
        // g[iC].x10[IdOUT]=((bne.pf1[IdBNE]+h12[j])+mp.df*(1-mv.pc)*mp.tpm12*(g[iC].v10[IdOUT]+mv.logsumK))/(1-mp.df*(1-mv.pc)*(1-mp.tpm12));
        // g[iC].x11[IdOUT]=g[iC].x10[IdOUT];
        
        // value functions for firm 2
        // f2=(mp.df*(1-mv.pc)*mp.tpm12)/(1-mp.df*(1-mv.pc)*(1-mp.tpm21));
        // g[iC].v20[IdOUT]=((bne.pf2[IdBNE]+h22[j])*(1+f2)+mp.df*(1-mv.pc)*(mp.tpm22+f2*mp.tpm21)*mv.logsumK)/(1-mp.df*(1-mv.pc)*(mp.tpm22+f2*mp.tpm21));
        // g[iC].v21[IdOUT]=g[iC].v20[IdOUT]-mv.kv;
        // g[iC].x20[IdOUT]=((bne.pf2[IdBNE]+h21[j])+mp.df*(1-mv.pc)*mp.tpm21*(g[iC].v20[IdOUT]+mv.logsumK))/(1-mp.df*(1-mv.pc)*(1-mp.tpm21));
        // g[iC].x21[IdOUT]=g[iC].x20[IdOUT];




        g[iC].p1[IdOUT]=f_logit(nF, mp.eta, NULL, g[iC].v11[IdOUT], g[iC].v10[IdOUT]);
        g[iC].p2[IdOUT]=f_logit(nF, mp.eta, NULL, g[iC].v21[IdOUT], g[iC].v20[IdOUT]);
        g[iC].seleqb[id(nC,ic1,ic2,iC,ieqb)]=1;
		
	    // expected costs: (c,c,c) corner game        
        //firm 1
        logsum1=f_lnsum(2, mp.eta, g[iC].v10[IdOUT], g[iC].v11[IdOUT]);		
        logsum2=f_lnsum(2, mp.eta, g[iC].v20[IdOUT], g[iC].v21[IdOUT]);		
		EP11 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v11[IdOUT])+mv.kv;
		EP10 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v10[IdOUT]);
		EP21 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v21[IdOUT])+mv.kv;
		EP20 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v20[IdOUT]);
        EPx1 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1];
		EPx2 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2];
        // -----------------------------------------------
        A1=   mp.tpm11*hc11[j] + mp.tpm21*hc12[j];
        B1=  (mp.df*(1-mv.pc))*mp.tpm11;
        C1=  (mp.df*(1-mv.pc))*mp.tpm21;
        A2=   mp.tpm12*hc11[j] + mp.tpm22*hc12[j];
        B2=  (mp.df*(1-mv.pc))*mp.tpm12;
        C2=  (mp.df*(1-mv.pc))*mp.tpm22;
        P1=  g[iC].p1[IdOUT];
        P2=  g[iC].p2[IdOUT];
        // CI_1(c,c,c,2): discounted cost of investing for firm 1, when it is firm 2's turn to invest (that is when firm 2 is investing)
        g[iC].ecx11[IdOUT] = -(A2 + EPx1 + A1*B2 - A2*B1 + B2*EP10 - B1*EPx1 - B2*EP10*P1 + B2*EP11*P1)  
                     /(B1 + C2 - B1*C2 + B2*C1 - 1); // CI2_ccc
        g[iC].ecv10[IdOUT]=  -(A1 + EP10 - A1*C2 + A2*C1 - C2*EP10 + C1*EPx1 + A1*C2*P2 - A2*C1*P2 
                                  + C1*g[iC].ecx11[IdOUT]*P2 - B1*EP10*P1 + B1*EP11*P1 + C2*EP10*P2 - C1*EPx1*P2 
                                  + B1*C2*EP10*P1 - B2*C1*EP10*P1 - B1*C2*EP11*P1 + B2*C1*EP11*P1 - B1*C2*EP10*P1*P2 
                                  + B2*C1*EP10*P1*P2 + B1*C2*EP11*P1*P2 - B2*C1*EP11*P1*P2)
                             /(B1 + C2 - B1*C2 + B2*C1 - C2*P2 + B1*C2*P2 - B2*C1*P2 - 1);// CN1_ccc
        g[iC].ecx10[IdOUT]=g[iC].ecx11[IdOUT];  // CN2_ccc
        g[iC].ecv11[IdOUT]=  -(A1 + EP11 + B1*g[iC].ecv10[IdOUT] + C1*g[iC].ecx11[IdOUT] - B1*g[iC].ecv10[IdOUT]*P1)/(B1*P1 - 1); //CI1_ccc     

        //firm2
        A1=   mp.tpm22*hc22[j] + mp.tpm12*hc21[j];
        B1=  (mp.df*(1-mv.pc))*mp.tpm22;
        C1=  (mp.df*(1-mv.pc))*mp.tpm12;
        A2=   mp.tpm21*hc22[j] + mp.tpm11*hc21[j];
        B2=  (mp.df*(1-mv.pc))*mp.tpm21;
        C2=  (mp.df*(1-mv.pc))*mp.tpm11;

        g[iC].ecx21[IdOUT] = -(A2 + EPx2 + A1*B2 - A2*B1 + B2*EP20 - B1*EPx2 - B2*EP20*P2 + B2*EP21*P2)  
                     /(B1 + C2 - B1*C2 + B2*C1 - 1); // CI2_ccc
        g[iC].ecv20[IdOUT]=  -(A1 + EP20 - A1*C2 + A2*C1 - C2*EP20 + C1*EPx2 + A1*C2*P1 - A2*C1*P1 
                                  + C1*g[iC].ecx21[IdOUT]*P1 - B1*EP20*P2 + B1*EP21*P2 + C2*EP20*P1 - C1*EPx2*P1 
                                  + B1*C2*EP20*P2 - B2*C1*EP20*P2 - B1*C2*EP21*P2 + B2*C1*EP21*P2 - B1*C2*EP20*P2*P1 
                                  + B2*C1*EP20*P2*P1 + B1*C2*EP21*P2*P1 - B2*C1*EP21*P2*P1)
                             /(B1 + C2 - B1*C2 + B2*C1 - C2*P1 + B1*C2*P1 - B2*C1*P1 - 1);// CN1_ccc
        g[iC].ecx20[IdOUT]=g[iC].ecx21[IdOUT];  // CN2_ccc
        g[iC].ecv21[IdOUT]=  -(A1 + EP21 + B1*g[iC].ecv20[IdOUT] + C1*g[iC].ecx21[IdOUT] - B1*g[iC].ecv20[IdOUT]*P2)/(B1*P2 - 1); //CI1_ccc

        /**************************************************************************/
        // SECTION 3.1: solve the (c1,c,c) edge game
        /**************************************************************************/
        ic2=iC; ieqb=0; //only one equilibrium in the edge  endgame
        
        // solve the (c1,c,c) end game
        mv.logsum1=f_lnsum(2, mp.eta,g[iC].v10[id(nC, iC, iC, iC, 0)],g[iC].v11[id(nC, iC, iC, iC, 0)]);
        mv.logsum2=f_lnsum(2, mp.eta,g[iC].v20[id(nC, iC, iC, iC, 0)],g[iC].v21[id(nC, iC, iC, iC, 0)]);
        
        for (ic1=iC+1;ic1<nC;ic1++) {
            if (((ic1z==ic1) && (ic2z==iC)) || (iC>iCz) || (iCz==0)) {

                j=(ic1-iC)*(nC-iC)+(ic2-iC);    // index used for storage in h
                
                g[iC].ic1[IdOUT]=ic1;
                g[iC].ic2[IdOUT]=ic2;
                g[iC].c1[IdOUT]=mp.cgrid[ic1];
                g[iC].c2[IdOUT]=mp.cgrid[ic2];
                g[iC].pf1[IdOUT]=bne.pf1[IdBNE];
                g[iC].pf2[IdOUT]=bne.pf2[IdBNE];
                g[iC].ieqb[IdOUT]=ieqb;
                g[iC].neqb[IdOUT]=1;
                g[iC].eqbtype[IdOUT]=1;   // 0 = pure strategy equilibrium
                g[iC].seleqb[IdOUT]=1;
                
                // solve for v10(c,0,0,m) and v11(c,0,0,m) for m=1,2
                f1=((1-mv.pc)*mp.df*mp.tpm21)/(1-mp.df*(1-mv.pc)*(1-mp.tpm12));
                g[iC].v11[IdOUT]=bne.pf1[IdBNE]-mv.kv+mp.df*(1-mv.pc)*(mp.tpm11*mv.logsum1+(mp.tpm21)*g[iC].x11[id(nC, iC, iC, iC, 0)])+h11[ih(nC,iC,iC,iC)];
                
                A  = bne.pf1[IdBNE]*(1+f1)+h11[j]+f1*h12[j];
                B  = mp.df*(1-mv.pc)*(mp.tpm11+mp.tpm12*f1);
                C  = g[iC].v11[IdOUT];
                
                g[iC].v10[IdOUT]=f_solvevf(A, B, C, mp.eta, &mp);
                g[iC].x11[IdOUT]=(bne.pf1[IdBNE]+h12[j]+mp.df*(1-mv.pc)*mp.tpm12*f_lnsum(2, mp.eta,g[iC].v11[IdOUT],g[iC].v10[IdOUT]))/(1-mp.df*(1-mv.pc)*mp.tpm22);
                g[iC].x10[IdOUT]=g[iC].x11[IdOUT];
                g[iC].p1[IdOUT]=f_logit(nF, mp.eta, NULL, g[iC].v11[IdOUT], g[iC].v10[IdOUT]);
                
                //solve for v20(c,0,0,m) and v21(c,0,0,m) for m=1,2
                f1=mp.tpm12*mp.df*(1-mv.pc)*(1-g[iC].p1[IdOUT])/(1-mp.df*(1-mv.pc)*mp.tpm11*(1-g[iC].p1[IdOUT]));
                g[iC].x21[IdOUT]=bne.pf2[IdBNE]+h21[ih(nC,iC,iC,iC)]+mp.df*(1-mv.pc)*(mp.tpm21*mv.logsum2+mp.tpm11*g[iC].x21[id(nC, iC, iC, iC, 0)]);   // BUG HERE ;)
                g[iC].v20[IdOUT]=(bne.pf2[IdBNE]*(1+f1)+h22[j]+f1*h21[j]+mp.df*(1-mv.pc)*(mp.tpm22+mp.tpm21*f1)*mv.logsumK + g[iC].x21[IdOUT]*mp.df*(1-mv.pc)*g[iC].p1[IdOUT]*(mp.tpm12+mp.tpm11*f1))
                /(1-mp.df*(1-mv.pc)*mp.tpm22-mp.df*(1-mv.pc)*mp.tpm21*f1);
                g[iC].v21[IdOUT]=g[iC].v20[IdOUT]-mv.kv;
                g[iC].x20[IdOUT]=(bne.pf2[IdBNE]+h21[j]+mp.df*(1-mv.pc)*(mp.tpm21*(g[iC].v20[IdOUT]+mv.logsumK)+mp.tpm11*g[iC].p1[IdOUT]*g[iC].x21[IdOUT]))
                /(1-mp.df*(1-mv.pc)*mp.tpm11*(1-g[iC].p1[IdOUT])); // BUG HERE
                g[iC].p2[IdOUT]=f_logit(nF, mp.eta, NULL, g[iC].v21[IdOUT], g[iC].v20[IdOUT]);
                
                // expected costs: (c1,c,c) corner game
                //common for both firms
                logsum1=f_lnsum(2, mp.eta, g[iC].v10[IdOUT], g[iC].v11[IdOUT]);     
                logsum2=f_lnsum(2, mp.eta, g[iC].v20[IdOUT], g[iC].v21[IdOUT]);     
                EP11 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v11[IdOUT])+mv.kv;
                EP10 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v10[IdOUT]);
                EP21 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v21[IdOUT])+mv.kv;
                EP20 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v20[IdOUT]);
                EPx1 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1];
                EPx2 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2];

                //firm 1 cost recursion
                A1=   mp.tpm11*hc11[j] + mp.tpm21*hc12[j];
                B1=  (mp.df*(1-mv.pc))*mp.tpm11;
                C1=  (mp.df*(1-mv.pc))*mp.tpm21;
                A2=   mp.tpm12*hc11[j] + mp.tpm22*hc12[j];
                B2=  (mp.df*(1-mv.pc))*mp.tpm12;
                C2=  (mp.df*(1-mv.pc))*mp.tpm22;
                P1=  g[iC].p1[IdOUT];
                P2=  g[iC].p2[IdOUT];
                
                 g[iC].ecv11[IdOUT]=
                         EP11 + mp.tpm11*hc11[ih(nC,iC,iC,iC)] + mp.tpm21*hc12[ih(nC,iC,iC,iC)]
                       +B1*(     g[iC].p1[id(nC, iC, iC, iC, 0)] *g[iC].ecv11[id(nC, iC, iC, iC, 0)]
                             +(1-g[iC].p1[id(nC, iC, iC, iC, 0)])*g[iC].ecv10[id(nC, iC, iC, iC, 0)])
                       +C1*(     g[iC].p2[id(nC, iC, iC, iC, 0)] *g[iC].ecx11[id(nC, iC, iC, iC, 0)]
                             +(1-g[iC].p2[id(nC, iC, iC, iC, 0)])*g[iC].ecx10[id(nC, iC, iC, iC, 0)]);
                CI1 = g[iC].ecv11[IdOUT];
                
                g[iC].ecv10[IdOUT]= -(A1 + EP10 - A1*C2 + A2*C1 - C2*EP10 + C1*EPx1 + B1*CI1*P1 - B1*C2*CI1*P1 + B2*C1*CI1*P1)
                                    /(B1 + C2 - B1*C2 + B2*C1 - B1*P1 + B1*C2*P1 - B2*C1*P1 - 1); // CN1_c1cc 
                g[iC].ecx10[IdOUT]= -(A2 + EPx1 + A1*B2 - A2*B1 + B2*EP10 - B1*EPx1 - A1*B2*P1 + A2*B1*P1 + B2*CI1*P1 - B2*EP10*P1 + B1*EPx1*P1)
                                    /(B1 + C2 - B1*C2 + B2*C1 - B1*P1 + B1*C2*P1 - B2*C1*P1 - 1); // CN2_c1cc 
                g[iC].ecx11[IdOUT]= -(A2 + EPx1 + A1*B2 - A2*B1 + B2*EP10 - B1*EPx1 - A1*B2*P1 + A2*B1*P1 + B2*CI1*P1 - B2*EP10*P1 + B1*EPx1*P1)
                                    /(B1 + C2 - B1*C2 + B2*C1 - B1*P1 + B1*C2*P1 - B2*C1*P1 - 1); // CI2_c1cc 

                //firm 2 cost recursion c1,c,c game (similar to firm 1 in c,c2,c game)
                A1=   mp.tpm22*hc22[j] + mp.tpm12*hc21[j];
                B1=  (mp.df*(1-mv.pc))*mp.tpm22;
                C1=  (mp.df*(1-mv.pc))*mp.tpm12;
                A2=   mp.tpm21*hc22[j] + mp.tpm11*hc21[j];
                B2=  (mp.df*(1-mv.pc))*mp.tpm21;
                C2=  (mp.df*(1-mv.pc))*mp.tpm11;
                
                g[iC].ecx21[IdOUT]=
                        EPx2 + mp.tpm21*hc22[ih(nC,iC,iC,iC)] + mp.tpm11*hc21[ih(nC,iC,iC,iC)]
                      +B2*(     g[iC].p2[id(nC, iC, iC, iC, 0)] *g[iC].ecv21[id(nC, iC, iC, iC, 0)]
                            +(1-g[iC].p2[id(nC, iC, iC, iC, 0)])*g[iC].ecv20[id(nC, iC, iC, iC, 0)])
                      +C2*(     g[iC].p1[id(nC, iC, iC, iC, 0)] *g[iC].ecx21[id(nC, iC, iC, iC, 0)]
                            +(1-g[iC].p1[id(nC, iC, iC, iC, 0)])*g[iC].ecx20[id(nC, iC, iC, iC, 0)]);

                CI2=g[iC].ecx21[IdOUT];
                g[iC].ecv20[IdOUT]= -(A1 + EP20 - A1*C2 + A2*C1 - C2*EP20 + C1*EPx2 + A1*C2*P1 - A2*C1*P1 + C1*CI2*P1 - B1*EP20*P2 
                                         + B1*EP21*P2 + C2*EP20*P1 - C1*EPx2*P1 + B1*C2*EP20*P2 - B2*C1*EP20*P2 - B1*C2*EP21*P2 
                                         + B2*C1*EP21*P2 - B1*C2*EP20*P1*P2 + B2*C1*EP20*P1*P2 + B1*C2*EP21*P1*P2 - B2*C1*EP21*P1*P2)
                                    /(B1 + C2 - B1*C2 + B2*C1 - C2*P1 + B1*C2*P1 - B2*C1*P1 - 1);
                g[iC].ecv21[IdOUT]= -(A1 + EP21 - A1*C2 + A2*C1 + B1*EP20 - B1*EP21 - C2*EP21 + C1*EPx2 - B1*C2*EP20 + B2*C1*EP20 
                                         + B1*C2*EP21 - B2*C1*EP21 + A1*C2*P1 - A2*C1*P1 + C1*CI2*P1 - B1*EP20*P2 + B1*EP21*P2 
                                         + C2*EP21*P1 - C1*EPx2*P1 + B1*C2*EP20*P1 - B2*C1*EP20*P1 + B1*C2*EP20*P2 - B1*C2*EP21*P1 
                                         - B2*C1*EP20*P2 + B2*C1*EP21*P1 - B1*C2*EP21*P2 + B2*C1*EP21*P2 - B1*C2*EP20*P1*P2 + B2*C1*EP20*P1*P2 
                                         + B1*C2*EP21*P1*P2 - B2*C1*EP21*P1*P2)
                                    /(B1 + C2 - B1*C2 + B2*C1 - C2*P1 + B1*C2*P1 - B2*C1*P1 - 1);
                g[iC].ecx20[IdOUT]= -(A2 + EPx2 + A1*B2 - A2*B1 + B2*EP20 - B1*EPx2 + C2*CI2*P1 - B2*EP20*P2 
                                                + B2*EP21*P2 - B1*C2*CI2*P1 + B2*C1*CI2*P1)
                                    /(B1 + C2 - B1*C2 + B2*C1 - C2*P1 + B1*C2*P1 - B2*C1*P1 - 1);
            }
//         printf("ic1=%d, ic2=%d, iC=%d, p1=%f, p2=%f, v21=%f, v20=%f\n", ic1, ic2, iC,g[iC].p1[IdOUT],g[iC].p2[IdOUT],g[iC].v21[IdOUT], g[iC].v20[IdOUT]);
        }
        
        /**************************************************************************/
        // SECTION 3.2: solve the (c,c2,c) edge game
        /**************************************************************************/
        ic1=iC; ieqb=0; //only one equilibrium in the edge endgame
        for (ic2=iC+1;ic2<nC;ic2++) {
            j=(ic1-iC)*(nC-iC)+(ic2-iC);    // index used for storage in h
            if (((ic2z==ic2) && (ic1z==iC)) || (iC>iCz) || (iCz==0)) {
                
                g[iC].ic1[IdOUT]=ic1;
                g[iC].ic2[IdOUT]=ic2;
                g[iC].c1[IdOUT]=mp.cgrid[ic1];
                g[iC].c2[IdOUT]=mp.cgrid[ic2];
                g[iC].pf1[IdOUT]=bne.pf1[IdBNE];
                g[iC].pf2[IdOUT]=bne.pf2[IdBNE];
                g[iC].ieqb[IdOUT]=ieqb;
                g[iC].neqb[IdOUT]=1;
                g[iC].eqbtype[IdOUT]=0;   // 0 = pure strategy equilibrium
                g[iC].seleqb[IdOUT]=1;
                
                // solve for v20(0,c,0,m) and v21(0,c,0,m) for m=1,2
                f2=(mp.df*(1-mv.pc)*mp.tpm12)/(1-mp.df*(1-mv.pc)*(1-mp.tpm21));
                g[iC].v21[IdOUT]=bne.pf2[IdBNE]-mv.kv+mp.df*(1-mv.pc)*(mp.tpm22*mv.logsum2+(mp.tpm12)*g[iC].x21[id(nC, iC, iC, iC, 0)])+h22[ih(nC,iC,iC,iC)];
                A  = bne.pf2[IdBNE]*(1+f2)+h22[j]+f2*h21[j];
                B  = mp.df*(1-mv.pc)*(mp.tpm22+mp.tpm21*f2);
                C  = g[iC].v21[IdOUT];
                
                g[iC].v20[IdOUT]=f_solvevf(A, B, C, mp.eta, &mp);
                g[iC].x21[IdOUT]=(bne.pf2[IdBNE]+h21[j]+mp.df*(1-mv.pc)*mp.tpm21*f_lnsum(2, mp.eta,g[iC].v21[IdOUT],g[iC].v20[IdOUT]))/(1-mp.df*(1-mv.pc)*mp.tpm11);
                g[iC].x20[IdOUT]=g[iC].x21[IdOUT];
                g[iC].p2[IdOUT]=f_logit(nF, mp.eta, NULL, g[iC].v21[IdOUT], g[iC].v20[IdOUT]);
                
                // solve for v10(0,c,0,m) and v21(0,c,0,m) for m=1,2
                f2=mp.tpm21*mp.df*(1-mv.pc)*(1-g[iC].p2[IdOUT])/(1-mp.df*(1-mv.pc)*mp.tpm22*(1-g[iC].p2[IdOUT]));
                g[iC].x11[IdOUT]=bne.pf1[IdBNE]+h12[ih(nC,iC,iC,iC)]+mp.df*(1-mv.pc)*(mp.tpm12*mv.logsum1+mp.tpm22*g[iC].x11[id(nC, iC, iC, iC, 0)]);         // BUG HERE
                g[iC].v10[IdOUT]=(bne.pf1[IdBNE]*(1+f2)+h11[j]+f2*h12[j]+mp.df*(1-mv.pc)*(mp.tpm11+mp.tpm12*f2)*mv.logsumK + g[iC].x11[IdOUT]*mp.df*(1-mv.pc)*g[iC].p2[IdOUT]*(mp.tpm21+mp.tpm22*f2))
                /(1-mp.df*(1-mv.pc)*mp.tpm11-mp.df*(1-mv.pc)*mp.tpm12*f2);
                g[iC].v11[IdOUT]=g[iC].v10[IdOUT]-mv.kv;
                g[iC].x10[IdOUT]=(bne.pf1[IdBNE]+h12[j]+mp.df*(1-mv.pc)*(mp.tpm12*(g[iC].v10[IdOUT]+mv.logsumK)+mp.tpm22*g[iC].p2[IdOUT]*g[iC].x11[IdOUT])) /(1-mp.df*(1-mv.pc)*mp.tpm22*(1-g[iC].p2[IdOUT]));
                g[iC].p1[IdOUT]=f_logit(nF, mp.eta, NULL, g[iC].v11[IdOUT], g[iC].v10[IdOUT]);
                
                
                // expected costs: (c,c2,c) edge game
                //common
                logsum1=f_lnsum(2, mp.eta, g[iC].v10[IdOUT], g[iC].v11[IdOUT]);     
                logsum2=f_lnsum(2, mp.eta, g[iC].v20[IdOUT], g[iC].v21[IdOUT]);     
                EP11 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v11[IdOUT])+mv.kv;
                EP10 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v10[IdOUT]);
                EP21 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v21[IdOUT])+mv.kv;
                EP20 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v20[IdOUT]);
                EPx1 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1];
                EPx2 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2];
                
                //firm 1
                A1=   mp.tpm11*hc11[j] + mp.tpm21*hc12[j];
                B1=  (mp.df*(1-mv.pc))*mp.tpm11;
                C1=  (mp.df*(1-mv.pc))*mp.tpm21;
                A2=   mp.tpm12*hc11[j] + mp.tpm22*hc12[j];
                B2=  (mp.df*(1-mv.pc))*mp.tpm12;
                C2=  (mp.df*(1-mv.pc))*mp.tpm22;
                P1=  g[iC].p1[IdOUT];
                P2=  g[iC].p2[IdOUT];
                
                g[iC].ecx11[IdOUT]=
                        EPx1 + mp.tpm12*hc11[ih(nC,iC,iC,iC)] + mp.tpm22*hc12[ih(nC,iC,iC,iC)]
                      +B2*(     g[iC].p1[id(nC, iC, iC, iC, 0)] *g[iC].ecv11[id(nC, iC, iC, iC, 0)]
                            +(1-g[iC].p1[id(nC, iC, iC, iC, 0)])*g[iC].ecv10[id(nC, iC, iC, iC, 0)])
                      +C2*(     g[iC].p2[id(nC, iC, iC, iC, 0)] *g[iC].ecx11[id(nC, iC, iC, iC, 0)]
                            +(1-g[iC].p2[id(nC, iC, iC, iC, 0)])*g[iC].ecx10[id(nC, iC, iC, iC, 0)]);

                CI2=g[iC].ecx11[IdOUT];
                g[iC].ecv10[IdOUT]= -(A1 + EP10 - A1*C2 + A2*C1 - C2*EP10 + C1*EPx1 + A1*C2*P2 - A2*C1*P2 
                                         + C1*CI2*P2 - B1*EP10*P1 + B1*EP11*P1 + C2*EP10*P2 - C1*EPx1*P2 
                                         + B1*C2*EP10*P1 - B2*C1*EP10*P1 - B1*C2*EP11*P1 + B2*C1*EP11*P1 - B1*C2*EP10*P1*P2 
                                         + B2*C1*EP10*P1*P2 + B1*C2*EP11*P1*P2 - B2*C1*EP11*P1*P2)
                                    /(B1 + C2 - B1*C2 + B2*C1 - C2*P2 + B1*C2*P2 - B2*C1*P2 - 1); // CN1_cc2c
                g[iC].ecv11[IdOUT]= -(A1 + EP11 - A1*C2 + A2*C1 + B1*EP10 - B1*EP11 - C2*EP11 + C1*EPx1 - B1*C2*EP10 
                                         + B2*C1*EP10 + B1*C2*EP11 - B2*C1*EP11 + A1*C2*P2 - A2*C1*P2 + C1*CI2*P2 - B1*EP10*P1 
                                         + B1*EP11*P1 + C2*EP11*P2 - C1*EPx1*P2 + B1*C2*EP10*P1 - B2*C1*EP10*P1 + B1*C2*EP10*P2 
                                         - B1*C2*EP11*P1 - B2*C1*EP10*P2 + B2*C1*EP11*P1 - B1*C2*EP11*P2 + B2*C1*EP11*P2 
                                         - B1*C2*EP10*P1*P2 + B2*C1*EP10*P1*P2 + B1*C2*EP11*P1*P2 - B2*C1*EP11*P1*P2)
                                    /(B1 + C2 - B1*C2 + B2*C1 - C2*P2 + B1*C2*P2 - B2*C1*P2 - 1); // CI1_cc2c  
                g[iC].ecx10[IdOUT]= -(A2 + EPx1 + A1*B2 - A2*B1 + B2*EP10 - B1*EPx1 + C2*CI2*P2 - B2*EP10*P1 
                                         + B2*EP11*P1 - B1*C2*CI2*P2 + B2*C1*CI2*P2)
                                    /(B1 + C2 - B1*C2 + B2*C1 - C2*P2 + B1*C2*P2 - B2*C1*P2 - 1); // CN2_cc2c

                //firm 2 in c,c2,c edge game (similar to firm 1 in c1,c,c)
                A1=   mp.tpm22*hc22[j] + mp.tpm12*hc21[j];
                B1=  (mp.df*(1-mv.pc))*mp.tpm22;
                C1=  (mp.df*(1-mv.pc))*mp.tpm12;
                A2=   mp.tpm21*hc22[j] + mp.tpm11*hc21[j];
                B2=  (mp.df*(1-mv.pc))*mp.tpm21;
                C2=  (mp.df*(1-mv.pc))*mp.tpm11;

                g[iC].ecv21[IdOUT]=
                          EP21 + mp.tpm22*hc22[ih(nC,iC,iC,iC)] + mp.tpm12*hc21[ih(nC,iC,iC,iC)]
                       +B1*(     g[iC].p2[id(nC, iC, iC, iC, 0)] *g[iC].ecv21[id(nC, iC, iC, iC, 0)]
                             +(1-g[iC].p2[id(nC, iC, iC, iC, 0)])*g[iC].ecv20[id(nC, iC, iC, iC, 0)])
                       +C1*(     g[iC].p1[id(nC, iC, iC, iC, 0)] *g[iC].ecx21[id(nC, iC, iC, iC, 0)]
                             +(1-g[iC].p1[id(nC, iC, iC, iC, 0)])*g[iC].ecx20[id(nC, iC, iC, iC, 0)]);
                CI1 = g[iC].ecv21[IdOUT];

                g[iC].ecv20[IdOUT]= -(A1 + EP20 - A1*C2 + A2*C1 - C2*EP20 + C1*EPx2 + B1*CI1*P2 - B1*C2*CI1*P2 + B2*C1*CI1*P2)
                                    /(B1 + C2 - B1*C2 + B2*C1 - B1*P2 + B1*C2*P2 - B2*C1*P2 - 1);
                g[iC].ecx20[IdOUT]= -(A2 + EPx2 + A1*B2 - A2*B1 + B2*EP20 - B1*EPx2 - A1*B2*P2 + A2*B1*P2 + B2*CI1*P2 - B2*EP20*P2 + B1*EPx2*P2)
                                    /(B1 + C2 - B1*C2 + B2*C1 - B1*P2 + B1*C2*P2 - B2*C1*P2 - 1);
                g[iC].ecx21[IdOUT]= -(A2 + EPx2 + A1*B2 - A2*B1 + B2*EP20 - B1*EPx2 - A1*B2*P2 + A2*B1*P2 + B2*CI1*P2 - B2*EP20*P2 + B1*EPx2*P2)
                                    /(B1 + C2 - B1*C2 + B2*C1 - B1*P2 + B1*C2*P2 - B2*C1*P2 - 1);
            }
        }
                
        /**************************************************************************/
        // SECTION 3.3: solve the (c1,c2,c) interior game
        /**************************************************************************/
        
        //loop over the interier points
        for (ic1=iC+1;ic1<nC;ic1++) {
            for (ic2=iC+1;ic2<nC;ic2++) {
                if (((ic1z==ic1) && (ic2z==ic2)) || (iC>iCz) || (iCz==0)) {
//                 if (true) {
                    //precompute for speed and convenience
                    j=(ic1-iC)*(nC-iC)+(ic2-iC);    // index used for storage in h
                    i1=id(nC, ic1, iC, iC, 0); //c1 c c
                    i2=id(nC, iC, ic2, iC, 0); //c c2 c
                    
                    // populate brc sturcture
                    brc.pf[0]=bne.pf1[IdBNE];
                    brc.pf[1]=bne.pf2[IdBNE];
                    brc.h11=h11[j];
                    brc.h12=h12[j];
                    brc.h21=h21[j];
                    brc.h22=h22[j];
                    brc.betnpc=mp.df*(1-mv.pc);
                    
                    //precompute values of investing
                    brc.v1[0]=bne.pf1[IdBNE]+h11[ih(nC,iC,ic2,iC)]-mv.kv+mp.df*(1-mv.pc)*mp.tpm11*f_lnsum(2, mp.eta, g[iC].v11[i2],g[iC].v10[i2])+mp.df*(1-mv.pc)*mp.tpm21*(g[iC].p2[i2]*g[iC].x11[i2]+(1-g[iC].p2[i2])*g[iC].x10[i2]);
                    brc.x1[0]=bne.pf1[IdBNE]+h12[ih(nC,ic1,iC,iC)]      +mp.df*(1-mv.pc)*mp.tpm12*f_lnsum(2, mp.eta, g[iC].v11[i1],g[iC].v10[i1])+mp.df*(1-mv.pc)*mp.tpm22*(g[iC].p2[i1]*g[iC].x11[i1]+(1-g[iC].p2[i1])*g[iC].x10[i1]);
                    brc.v1[1]=bne.pf2[IdBNE]+h22[ih(nC,ic1,iC,iC)]-mv.kv+mp.df*(1-mv.pc)*mp.tpm22*f_lnsum(2, mp.eta, g[iC].v21[i1],g[iC].v20[i1])+mp.df*(1-mv.pc)*mp.tpm12*(g[iC].p1[i1]*g[iC].x21[i1]+(1-g[iC].p1[i1])*g[iC].x20[i1]);
                    brc.x1[1]=bne.pf2[IdBNE]+h21[ih(nC,iC,ic2,iC)]      +mp.df*(1-mv.pc)*mp.tpm21*f_lnsum(2, mp.eta, g[iC].v21[i2],g[iC].v20[i2])+mp.df*(1-mv.pc)*mp.tpm11*(g[iC].p1[i2]*g[iC].x21[i2]+(1-g[iC].p1[i2])*g[iC].x20[i2]);
                    
                    
                    //initialize the number of equilibria to be stored in the g structure
                    g[iC].neqb[id(nC,ic1,ic2,iC,0)]=0;
                    
                    //find equilibria
                    if (mp.eta>0){
                        s_ffxp((f_br2p1_am), &mp, &brc, &g[iC], ic1, ic2, iC, nC);
                    }
                    else {
                        analytical_am(&mp, &brc, &g[iC], ic1, ic2, iC, nC);
                    }
                    
                    //number of equilibria found
                    //neqb=f_neqb(nC, ic1, ic2, iC, g); XXXX
                    neqb=(int)g[iC].neqb[id(nC,ic1,ic2,iC,0)];
                    
                    // Save indices and costs (remaining elements of g are stored inside f_br2p1_am
                    for (ieqb=0;ieqb<neqb;ieqb++) {
                        
                        g[iC].ic1[IdOUT]=ic1;               //0
                        g[iC].ic2[IdOUT]=ic2;               //1
                        g[iC].ieqb[IdOUT]=ieqb;             //2
                        g[iC].c1[IdOUT]=mp.cgrid[ic1];      //3
                        g[iC].c2[IdOUT]=mp.cgrid[ic2];      //4
                        // 5-11: eqbtype, v10,v11,b20, v21, p1, p2: stored inside f_br2p1_am
//                         if ((mp.eta==0) && fabs(mp.tpm11)<mp.ctol && fabs(mp.tpm22)<mp.ctol) {
//                             g[iC].seleqb[id(nC,ic1,ic2,iC,ieqb)]=(ieqb==0);  //12
//                         }
//                         else {
                            g[iC].seleqb[id(nC,ic1,ic2,iC,ieqb)]=(ieqb==f_ers(&mp, iC, ic1, ic2, neqb,eqstring)?1.0:0.0);  //12
//                         }
                        g[iC].pf1[IdOUT]=bne.pf1[IdBNE];           //13
                        g[iC].pf2[IdOUT]=bne.pf2[IdBNE];           //14
                                
                        // expected costs: (c1,c2,c) corner game
                        //firm 1
                        logsum1=f_lnsum(2, mp.eta, g[iC].v10[IdOUT], g[iC].v11[IdOUT]);     
                        logsum2=f_lnsum(2, mp.eta, g[iC].v20[IdOUT], g[iC].v21[IdOUT]);     
                        EP11 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v11[IdOUT])+mv.kv;
                        EP10 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1]-mp.eta*(EULER+logsum1-g[iC].v10[IdOUT]);
                        EP21 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v21[IdOUT])+mv.kv;
                        EP20 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2]-mp.eta*(EULER+logsum2-g[iC].v20[IdOUT]);
                        EPx1 = mp.dt*bne.s1[IdBNE]*mp.cgrid[ic1];
                        EPx2 = mp.dt*bne.s2[IdBNE]*mp.cgrid[ic2];

                                    
                        g[iC].ecx11[IdOUT]=
                                EPx1 + mp.tpm12*hc11[ih(nC,ic1,iC,iC)] + mp.tpm22*hc12[ih(nC,ic1,iC,iC)]
                              +(mp.df*(1-mv.pc))*mp.tpm12*(
                                        g[iC].p1[id(nC, ic1, iC, iC, 0)] *g[iC].ecv11[id(nC, ic1, iC, iC, 0)]
                                    +(1-g[iC].p1[id(nC, ic1, iC, iC, 0)])*g[iC].ecv10[id(nC, ic1, iC, iC, 0)])
                              +(mp.df*(1-mv.pc))*mp.tpm22*(
                                        g[iC].p2[id(nC, ic1, iC, iC, 0)] *g[iC].ecx11[id(nC, ic1, iC, iC, 0)]
                                    +(1-g[iC].p2[id(nC, ic1, iC, iC, 0)])*g[iC].ecx10[id(nC, ic1, iC, iC, 0)]);
                
                        g[iC].ecv11[IdOUT]= 
                                  EP11 + mp.tpm11*hc11[ih(nC,iC,ic2,iC)] + mp.tpm21*hc12[ih(nC,iC,ic2,iC)] 
                              +(mp.df*(1-mv.pc))*mp.tpm11*(
                                        g[iC].p1[id(nC, iC, ic2, iC, 0)] *g[iC].ecv11[id(nC, iC, ic2, iC, 0)]
                                    +(1-g[iC].p1[id(nC, iC, ic2, iC, 0)])*g[iC].ecv10[id(nC, iC, ic2, iC, 0)])
                              +(mp.df*(1-mv.pc))*mp.tpm21*(
                                        g[iC].p2[id(nC, iC, ic2, iC, 0)] *g[iC].ecx11[id(nC, iC, ic2, iC, 0)]
                                    +(1-g[iC].p2[id(nC, iC, ic2, iC, 0)])*g[iC].ecx10[id(nC, iC, ic2, iC, 0)]);
                    
                        A1=   EP10 + mp.tpm11*hc11[j] + mp.tpm21*hc12[j]
                                   +(mp.df*(1-mv.pc))*mp.tpm11*g[iC].p1[IdOUT]*g[iC].ecv11[IdOUT]+
                                   +(mp.df*(1-mv.pc))*mp.tpm21*g[iC].p2[IdOUT]*g[iC].ecx11[IdOUT];
                        B1=  (mp.df*(1-mv.pc))*mp.tpm11*(1-g[iC].p1[IdOUT]);
                        C1=  (mp.df*(1-mv.pc))*mp.tpm21*(1-g[iC].p2[IdOUT]);
                        A2=  EPx1 + mp.tpm12*hc11[j] + mp.tpm22*hc12[j]
                                  +(mp.df*(1-mv.pc))*mp.tpm12*g[iC].p1[IdOUT]*g[iC].ecv11[IdOUT]+
                                  +(mp.df*(1-mv.pc))*mp.tpm22*g[iC].p2[IdOUT]*g[iC].ecx11[IdOUT];
                        B2=  (mp.df*(1-mv.pc))*mp.tpm12*(1-g[iC].p1[IdOUT]);
                        C2=  (mp.df*(1-mv.pc))*mp.tpm22*(1-g[iC].p2[IdOUT]);

                        g[iC].ecv10[IdOUT] =   (A1*(1-C2)+C1*A2) 
                                             / ((1-B1)*(1-C2)-C1*B2);
                        g[iC].ecx10[IdOUT] =   (A2+B2*g[iC].ecv10[IdOUT])
                                             / (1-C2);

                         //firm 2

                        g[iC].ecx21[IdOUT]=
                                EPx2 + mp.tpm21*hc22[ih(nC,iC,ic2,iC)] + mp.tpm11*hc21[ih(nC,iC,ic2,iC)]
                              +(mp.df*(1-mv.pc))*mp.tpm21*(
                                        g[iC].p2[id(nC, iC, ic2, iC, 0)] *g[iC].ecv21[id(nC, iC, ic2, iC, 0)]
                                    +(1-g[iC].p2[id(nC, iC, ic2, iC, 0)])*g[iC].ecv20[id(nC, iC, ic2, iC, 0)])
                              +(mp.df*(1-mv.pc))*mp.tpm11*(
                                        g[iC].p1[id(nC, iC, ic2, iC, 0)] *g[iC].ecx21[id(nC, iC, ic2, iC, 0)]
                                    +(1-g[iC].p1[id(nC, iC, ic2, iC, 0)])*g[iC].ecx20[id(nC, iC, ic2, iC, 0)]);
                        g[iC].ecv21[IdOUT]= 
                                  EP21 + mp.tpm22*hc22[ih(nC,ic1,iC,iC)] + mp.tpm12*hc21[ih(nC,ic1,iC,iC)] 
                              +(mp.df*(1-mv.pc))*mp.tpm22*(
                                        g[iC].p2[id(nC, ic1, iC, iC, 0)] *g[iC].ecv21[id(nC, ic1, iC, iC, 0)]
                                    +(1-g[iC].p2[id(nC, ic1, iC, iC, 0)])*g[iC].ecv20[id(nC, ic1, iC, iC, 0)])
                              +(mp.df*(1-mv.pc))*mp.tpm12*(
                                        g[iC].p1[id(nC, ic1, iC, iC, 0)] *g[iC].ecx21[id(nC, ic1, iC, iC, 0)]
                                    +(1-g[iC].p1[id(nC, ic1, iC, iC, 0)])*g[iC].ecx20[id(nC, ic1, iC, iC, 0)]);

                        A1=   EP20 + mp.tpm22*hc22[j] + mp.tpm12*hc21[j]
                                   +(mp.df*(1-mv.pc))*mp.tpm22*g[iC].p2[IdOUT]*g[iC].ecv21[IdOUT]+
                                   +(mp.df*(1-mv.pc))*mp.tpm12*g[iC].p1[IdOUT]*g[iC].ecx21[IdOUT];
                        B1=  (mp.df*(1-mv.pc))*mp.tpm22*(1-g[iC].p2[IdOUT]);
                        C1=  (mp.df*(1-mv.pc))*mp.tpm12*(1-g[iC].p1[IdOUT]);
                        A2=  EPx2 + mp.tpm21*hc22[j] + mp.tpm11*hc21[j]
                                  +(mp.df*(1-mv.pc))*mp.tpm21*g[iC].p2[IdOUT]*g[iC].ecv21[IdOUT]+
                                  +(mp.df*(1-mv.pc))*mp.tpm11*g[iC].p1[IdOUT]*g[iC].ecx21[IdOUT];
                        B2=  (mp.df*(1-mv.pc))*mp.tpm21*(1-g[iC].p2[IdOUT]);
                        C2=  (mp.df*(1-mv.pc))*mp.tpm11*(1-g[iC].p1[IdOUT]);
                    
                        g[iC].ecv20[IdOUT] =   (A1*(1-C2)+C1*A2) 
                                             / ((1-B1)*(1-C2)-C1*B2);
                        g[iC].ecx20[IdOUT] =   (A2+B2*g[iC].ecv20[IdOUT])
                                             / (1-C2);
                    }
                }
            }
        }

        //compute totals for the whole layer
        ieqb=0; //only one equilibrium in the edge endgame
        for (ic1=iC;ic1<nC;ic1++) {
            for (ic2=iC;ic2<nC;ic2++) {
                if (((ic1z==ic1) && (ic2z==ic2)) || (iC>iCz) || (iCz==0)) {
                    if (ic1>iC && ic2>iC) neqb=(int)g[iC].neqb[id(nC,ic1,ic2,iC,0)]; //interior
                    else neqb=1; //corner and edges 
                    for (ieqb=0;ieqb<neqb;ieqb++) {
                        // expected total costs from both firms (conditional on who has the right of move)
                        //g[iC].ecm1[IdOUT]=g[iC].p1[IdOUT]*(g[iC].ecv11[IdOUT]+g[iC].ecx21[IdOUT]) + (1-g[iC].p1[IdOUT])*(g[iC].ecv10[IdOUT]+g[iC].ecx20[IdOUT]);
                        //g[iC].ecm2[IdOUT]=g[iC].p2[IdOUT]*(g[iC].ecx11[IdOUT]+g[iC].ecv21[IdOUT]) + (1-g[iC].p2[IdOUT])*(g[iC].ecx10[IdOUT]+g[iC].ecv20[IdOUT]);
                        //social surplus = total revenue - total cost (conditional on who has the right of move)
                        g[iC].ecm1[IdOUT]=mp.dt*MAX(mp.cgrid[ic1],mp.cgrid[ic2])/(1-mp.df) - (g[iC].p1[IdOUT]*(g[iC].ecv11[IdOUT]+g[iC].ecx21[IdOUT]) + (1-g[iC].p1[IdOUT])*(g[iC].ecv10[IdOUT]+g[iC].ecx20[IdOUT]));
                        g[iC].ecm2[IdOUT]=mp.dt*MAX(mp.cgrid[ic1],mp.cgrid[ic2])/(1-mp.df) - (g[iC].p2[IdOUT]*(g[iC].ecx11[IdOUT]+g[iC].ecv21[IdOUT]) + (1-g[iC].p2[IdOUT])*(g[iC].ecx10[IdOUT]+g[iC].ecv20[IdOUT]));
                        // If run for monopoly, this is efficiency measure (revenue - cost)/profits
                        // g[iC].ecm1[IdOUT]=( mp.dt*MAX(mp.cgrid[ic1],mp.cgrid[ic2])/(1-mp.df) - (g[iC].p1[IdOUT]*(g[iC].ecv11[IdOUT]+g[iC].ecx21[IdOUT]) + (1-g[iC].p1[IdOUT])*(g[iC].ecv10[IdOUT]+g[iC].ecx20[IdOUT])) )
                        //                   /(g[iC].p1[IdOUT]*(g[iC].v11[IdOUT]) + (1-g[iC].p1[IdOUT])*(g[iC].v10[IdOUT]));
                        // g[iC].ecm2[IdOUT]=( mp.dt*MAX(mp.cgrid[ic1],mp.cgrid[ic2])/(1-mp.df) - (g[iC].p1[IdOUT]*(g[iC].ecv11[IdOUT]+g[iC].ecx21[IdOUT]) + (1-g[iC].p1[IdOUT])*(g[iC].ecv10[IdOUT]+g[iC].ecx20[IdOUT])) )
                        //                   /(g[iC].p1[IdOUT]*(g[iC].v11[IdOUT]+g[iC].x21[IdOUT]) + (1-g[iC].p1[IdOUT])*(g[iC].v10[IdOUT]+g[iC].x20[IdOUT]));
                    }
                }
            }
        }

//     }  // end parallel execution
    } //continue loop over iC
    
    free(h11);
    free(h21);
    free(h12);
    free(h22);

    free(hc11);
    free(hc21);
    free(hc12);
    free(hc22);
free(mv.eqbinfo);
}

double f_brp1_am(double p2, MPstruct* mp, Brc_amstruct* brc) {
    brc[0].f[0]=brc[0].betnpc*mp[0].tpm21*(1-p2)/(1-brc[0].betnpc*mp[0].tpm22*(1-p2));
    brc[0].A[0] = (1+brc[0].f[0])*brc[0].pf[0] +brc[0].betnpc*(mp[0].tpm21+mp[0].tpm22*brc[0].f[0])*p2*brc[0].x1[0]
            +brc[0].h11+brc[0].f[0]*brc[0].h12;
    brc[0].B[0] = brc[0].betnpc*(mp[0].tpm11+mp[0].tpm12*brc[0].f[0]);
    brc[0].v0[0]=f_solvevf(brc[0].A[0], brc[0].B[0], brc[0].v1[0], mp[0].eta, mp);
    brc[0].p[0]=f_logit(nF, mp[0].eta, NULL, brc[0].v1[0], brc[0].v0[0]);
    return brc[0].p[0];
}

double f_brp2_am(double p1, MPstruct* mp, Brc_amstruct* brc) {
    brc[0].f[1]=brc[0].betnpc*mp[0].tpm12*(1-p1)/(1-brc[0].betnpc*mp[0].tpm11*(1-p1));
    brc[0].A[1] = (1+brc[0].f[1])*brc[0].pf[1]+ brc[0].betnpc*(mp[0].tpm12+mp[0].tpm11*brc[0].f[1])*p1*brc[0].x1[1]
            +brc[0].h22+brc[0].f[1]*brc[0].h21;
    brc[0].B[1] = brc[0].betnpc*(mp[0].tpm22+mp[0].tpm21*brc[0].f[1]);
    brc[0].v0[1]=f_solvevf(brc[0].A[1], brc[0].B[1], brc[0].v1[1], mp[0].eta, mp);
    brc[0].p[1]=f_logit(nF, mp[0].eta, NULL, brc[0].v1[1], brc[0].v0[1]);
    return brc[0].p[1];
}

double f_br2p1_am(double p1, int ieqb, int eqbtype, MPstruct* mp, Brc_amstruct* brc, Gamestruct* gc, int ic1, int ic2, int iC, int nC) {
    double p2p1;
    int idz;
    p2p1=f_brp2_am(p1, mp, brc);
    p2p1=f_brp1_am(p2p1, mp, brc);
    if (ieqb>=0){
        idz=IdOUT;
        gc[0].ieqb[idz]=ieqb;          // already populated
        gc[0].neqb[id(nC,ic1,ic2,iC,0)]++;
        gc[0].eqbtype[idz]=eqbtype;
        gc[0].v10[idz]=brc[0].v0[0];
        gc[0].v11[idz]=brc[0].v1[0];
        gc[0].v20[idz]=brc[0].v0[1];
        gc[0].v21[idz]=brc[0].v1[1];
        gc[0].p1[idz]=brc[0].p[0];
        gc[0].p2[idz]=brc[0].p[1];
        gc[0].x11[idz]=brc[0].x1[0];
        gc[0].x21[idz]=brc[0].x1[1];
        gc[0].x10[idz]=(brc[0].pf[0]+brc[0].betnpc*mp[0].tpm12*f_lnsum(2, mp[0].eta, brc[0].v1[0],brc[0].v0[0])+brc[0].betnpc*mp[0].tpm22*brc[0].p[1]*brc[0].x1[0]+brc[0].h12)
        /(1-brc[0].betnpc*mp[0].tpm22*(1-brc[0].p[1]));
        gc[0].x20[idz]=(brc[0].pf[1]+brc[0].betnpc*mp[0].tpm21*f_lnsum(2, mp[0].eta, brc[0].v1[1],brc[0].v0[1])+brc[0].betnpc*mp[0].tpm11*brc[0].p[0]*brc[0].x1[1]+brc[0].h21)
        /(1-brc[0].betnpc*mp[0].tpm11*(1-brc[0].p[0]));
    }
    return p2p1;
}

void s_plot_bpr(MPstruct* mp, Brc_amstruct* brc, int ic1, int ic2, int iC) {
    mxArray *Mvars[4];
    const int Mpoints=1000;
    double *Mvar;
    char Mstring[500];
    int i, j;
    
    /*Plot best responce function*/
    for (i=0;i<4-1;i++) {
        Mvars[i]=mxCreateDoubleMatrix(Mpoints,1,mxREAL);
        Mvar=mxGetPr(Mvars[i]);
        for (j=0;j<Mpoints;Mvar[j++]=(double)j/Mpoints);
    }
    Mvars[3]=mxCreateDoubleMatrix(Mpoints,1,mxREAL);
    Mvar=mxGetPr(Mvars[3]);
    for (j=0;j<Mpoints;Mvar[j++]=f_brp1_am((double)j/Mpoints, mp, brc));
    mexCallMATLAB(0,NULL,4,Mvars,"plot");
    mexEvalString("hold on");
    
    Mvars[0]=mxCreateDoubleMatrix(Mpoints,1,mxREAL);
    Mvar=mxGetPr(Mvars[0]);
    for (j=0;j<Mpoints;Mvar[j++]=f_brp2_am((double)j/Mpoints, mp, brc));
    
    Mvars[3]=mxCreateDoubleMatrix(Mpoints,1,mxREAL);
    Mvar=mxGetPr(Mvars[3]);
    for (j=0;j<Mpoints;Mvar[j++]=(double)j/Mpoints);
    
    mexCallMATLAB(0,NULL,4,Mvars,"plot");
    mexEvalString("hold on");
    
    sprintf(Mstring,"title('c1=%1.1f c2=%1.1f c=%1.1f')",mp[0].cgrid[ic1],mp[0].cgrid[ic2], mp[0].cgrid[iC]);
    mexEvalString(Mstring);
    mexEvalString("xlabel('p1');");
    mexEvalString("ylabel('p2');");
    mexEvalString("hold on");
// if (ic2==nC-1) {
    mexEvalString("drawnow");
    mexEvalString("pause");
//     mexEvalString("input ('Press Enter to continue..');");
//     mexEvalString("close gcf;");
// }
}


