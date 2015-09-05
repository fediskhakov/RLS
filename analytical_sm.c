typedef struct brcstruct{
        double a0[2];
        double a1[2];
        double a2[2];
        int iF;
        double eta;
        double B1;
        double r1;
        double r2;
    } Brcstruct;  //shorthand for the arguments of the functions below    

void s_domainbr2(Brcstruct *brc,MPstruct *mp,double* domain);
double s_findeqb_bracketing(Brcstruct *brc,MPstruct *mp,double* domain);
int sa(double* x1, double x0, double l, double r,Brcstruct* brc,MPstruct *mp);
void findeqb_sa(Brcstruct* brc,MPstruct *mp,MVstruct *mv,double* domain);
double f_invbr(double p_i, Brcstruct* brc);
int f_invbr2(double* invbr2, double p_i, Brcstruct* brc);


static int VERBOSE;
void analytical_solutions (int ic1,int ic2,int iC,MPstruct mp,MVstruct mv,Bnestruct bne,Gamestruct* g,Brstruct* br, int* eqstring) {
    /*--------------------------------------------------------------------------------------
     * This function calculates all equilibria for given c1 c2 c
     * and puts them in proper output structures gc and brc
     *--------------------------------------------------------------------------------------*/
    int i, j, i1, i2, iF, nC;
    int ieqb, neqb;
    int save;
    double a0[2], a1[2], a2[2], pstar[4][2], A[2], B[2], v0_0[2], v0_1[2], pbr[2]; 
    int ir;
    double r1[4]={1,-1,-1,1}, r2[4]={1,1,-1,-1};
    double domain[2]; //domain of the second order best responce function
    Brcstruct brc; //shorthand for parameters
    double fxpoint;
    int ieqb_offset;
        int nPstar[2];

    nC=mp.nC;
    
    
    // ********************************************
    // Coefficients
    // ********************************************
    //precompute for speed and convenience
    i1=id(nC, ic1, iC, iC, 0); //c1 c c
    i2=id(nC, iC, ic2, iC, 0); //c c2 c
    
    //start filling out br output
    br[iC].ic1[brOUT]=ic1;
    br[iC].ic2[brOUT]=ic2;
    br[iC].c1[brOUT]=mp.cgrid[ic1];
    br[iC].c2[brOUT]=mp.cgrid[ic2];
    
    //Coefficients to A    
    br[iC].A0[0][brOUT]=bne.pf1[IdBNE]+mp.df*mv.pc*mv.h1[ih(nC, ic1, ic2, iC)];
    br[iC].A0[1][brOUT]=bne.pf2[IdBNE]+mp.df*mv.pc*mv.h2[ih(nC, ic1, ic2, iC)];

    br[iC].A1[0][brOUT]=mp.df*(
            mv.pc*mv.h1[ih(nC, ic1, iC, iC)]  + (1-mv.pc)*f_lnsum(nF, mp.eta, g[iC].v10[i1], g[iC].v11[i1])
            -mv.pc*mv.h1[ih(nC, ic1, ic2, iC)]
            );
    br[iC].A1[1][brOUT]=mp.df*(
            mv.pc*mv.h2[ih(nC, iC, ic2, iC)]  + (1-mv.pc)*f_lnsum(nF, mp.eta, g[iC].v20[i2], g[iC].v21[i2])
            -mv.pc*mv.h2[ih(nC, ic1, ic2, iC)]
            );

    //Coefficients to B
    br[iC].B1[brOUT]=mp.df*(1-mv.pc);

    //Coefficients to value of investing
    br[iC].C0[0][brOUT]=bne.pf1[IdBNE]-mv.kv
            + mp.df*( mv.pc*mv.h1[ih(nC, iC, ic2, iC)] + (1-mv.pc)*f_lnsum(nF, mp.eta, g[iC].v10[i2],g[iC].v11[i2])) ;
    br[iC].C0[1][brOUT]=bne.pf2[IdBNE]-mv.kv
            + mp.df*( mv.pc*mv.h2[ih(nC, ic1, iC, iC)] + (1-mv.pc)*f_lnsum(nF, mp.eta, g[iC].v20[i1],g[iC].v21[i1])) ;
    br[iC].C1[0][brOUT]=mp.df*(
             ( mv.pc*mv.h1[ih(nC, iC, iC, iC)]  + (1-mv.pc)*f_lnsum(nF, mp.eta, g[iC].v10[id(nC, iC, iC, iC, 0)],g[iC].v11[id(nC, iC, iC, iC, 0)])  ) -
             ( mv.pc*mv.h1[ih(nC, iC, ic2, iC)] + (1-mv.pc)*f_lnsum(nF, mp.eta, g[iC].v10[i2],g[iC].v11[i2])  )
            );
    br[iC].C1[1][brOUT]=mp.df*(
             ( mv.pc*mv.h2[ih(nC, iC, iC, iC)]  + (1-mv.pc)*f_lnsum(nF, mp.eta, g[iC].v20[id(nC, iC, iC, iC, 0)],g[iC].v21[id(nC, iC, iC, iC, 0)])  ) -
             ( mv.pc*mv.h2[ih(nC, ic1, iC, iC)] + (1-mv.pc)*f_lnsum(nF, mp.eta, g[iC].v20[i1],g[iC].v21[i1])  )
             );
   
    //Coefficients for polinomials
    for (iF=0;iF<2;iF++)
    {
        br[iC].D0[iF][brOUT]=br[iC].A0[iF][brOUT]-br[iC].C0[iF][brOUT]+br[iC].C0[iF][brOUT]*br[iC].B1[brOUT];
        br[iC].D1[iF][brOUT]=(br[iC].A1[iF][brOUT]-br[iC].C0[iF][brOUT]*br[iC].B1[brOUT]+(br[iC].B1[brOUT]-1)*br[iC].C1[iF][brOUT]);
        br[iC].D2[iF][brOUT]=-br[iC].B1[brOUT]*br[iC].C1[iF][brOUT];
    }


    // ********************************************
    // Switch by eta
    // ********************************************
    if (mp.eta==0.0)
    {   //Simple case when polinomials have constant coefficients
//         for (iF=0;iF<2;iF++) {
//             // Coefficients in polynomials short notation (Fedor agreed to this at 15:08 1.12.2011)
//             a0[iF]=br[iC].D0[iF][brOUT];
//             a1[iF]=br[iC].D1[iF][brOUT];
//             a2[iF]=br[iC].D2[iF][brOUT];
//             // Equilibrium candidates, NOTE iF index is the opponent here !!!!!!
//             pstar[0][iF]=0;
//             pstar[1][iF]=-(1/(2*a2[iF]))*(a1[iF]+sqrt(a1[iF]*a1[iF]-4*a0[iF]*a2[iF]));//positive root
//             pstar[2][iF]=-(1/(2*a2[iF]))*(a1[iF]-sqrt(a1[iF]*a1[iF]-4*a0[iF]*a2[iF]));//negative root
//             pstar[3][iF]=1;
// printf("iF=%d a0=%f a1=%f a2=%f roots= %f %f \n",iF,a0[iF],a1[iF],a2[iF],pstar[1][iF],pstar[2][iF]);
//         };
        
        for (iF=0;iF<2;iF++) {
            // Coefficients in polynomials short notation (Fedor agreed to this at 15:08 1.12.2011)
            a0[iF]=br[iC].D0[iF][brOUT];
            a1[iF]=br[iC].D1[iF][brOUT];
            a2[iF]=br[iC].D2[iF][brOUT];
        };
        for (iF=0;iF<2;iF++) {
            if (fabs(a2[iF])<mp.ctol) {
                pstar[0][iF]=0;
                pstar[1][iF]=-a0[iF]/a1[iF]; // single root
                pstar[2][iF]=1;
                nPstar[iF]=3;
                //printf("iF=%d a0=%f a1=%f a2=%f root SINGLE= %f 1.0\n",iF,a0[iF],a1[iF],a2[iF],pstar[1][iF]);                
            }
            else {
                // Equilibrium candidates, NOTE iF index is the opponent herD   !!!!!!
                pstar[0][iF]=0;
                pstar[1][iF]=-(1/(2*a2[iF]))*(a1[iF]+sqrt(a1[iF]*a1[iF]-4*a0[iF]*a2[iF]));//positive root
                pstar[2][iF]=-(1/(2*a2[iF]))*(a1[iF]-sqrt(a1[iF]*a1[iF]-4*a0[iF]*a2[iF]));//negative root
                pstar[3][iF]=1;
                nPstar[iF]=4;
                //printf("iF=%d a0=%f a1=%f a2=%f roots= %f %f \n",iF,a0[iF],a1[iF],a2[iF],pstar[1][iF],pstar[2][iF]);
            }
                        // Equilibrium candidates, NOTE iF index is the opponent here !!!!!!
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
                    pbr[0]=((a0[0]+a1[0]*pstar[i][0]+a2[0]*pstar[i][0]*pstar[i][0])<0);
                    //best responce for firm 2:
                    pbr[1]=((a0[1]+a1[1]*pstar[j][1]+a2[1]*pstar[j][1]*pstar[j][1])<0);
                    //check if probs really match
                    
                    if (fabs(pbr[0]-pstar[j][1])<mp.ctol && fabs(pbr[1]-pstar[i][0])<mp.ctol){
                        ieqb=ieqb+1;
                        save=1;
                    }
                }
                if (((i==1) || (i==(nPstar[0]-2))) && ((j==1) || j==(nPstar[1]-2))) { // mixed strategy
//                         printf("pstar[i=%d][0]=%f, pstar[j=%d][1]=%f\n", i,pstar[i][0],j,pstar[i][0]);
                    
                    if (((pstar[i][0]>=0) && (pstar[i][0]<=1)) && ((pstar[j][1]>=0) && (pstar[j][1]<=1))) {
                        ieqb=ieqb+1;
                        save=2;
                    }
                }
                
//                                 
//                 printf("iC=%d ic1=%d ic2=%d\nsave=%d ctol=%e pstar= ",iC,ic1,ic2,save,mp.ctol);
//                     printf("%f %f ", pstar[j][0],pstar[j][1]);
//                 printf("\nfabs(pbr[0]-pstar[j][1]=%f,   fabs(pbr[1]-pstar[i][0]))=%f", fabs(pbr[0]-pstar[j][1]),   fabs(pbr[1]-pstar[i][0]));
//                 printf("\npstar[i][0]=%f,  pstar[j][1]=%f", pstar[i][0],  pstar[j][1]);
//                 printf("\na0[0]=%f, a1[0]=%f,  a2[0]=%f", a0[0], a1[0],  a2[0]);
//                 printf("\na0[1]=%f, a1[1]=%f,  a2[1]=%f", a0[1], a1[1],  a2[1]);
//                 printf("\n\n");
// 
                
                if (ieqb>=MAXEQB) { 
                   printf("iC=%d ic1=%d ic2=%d\nsave=%d pstar=",iC,ic1,ic2,save);
                    for (j=0;j<4;j++) {
                     printf("%f %f ", pstar[j][0],pstar[j][1]);  
                    }
                   mexErrMsgTxt("Error in analytical.c, eta==0: number of equilibria exceeded MAXEQB");
                }
                //save new equilibrium (note IsOUT depends on ieqb)
                if (save>0) {
                    g[iC].ic1[IdOUT]=ic1;
                    g[iC].ic2[IdOUT]=ic2;
                    g[iC].pf1[IdOUT]=bne.pf1[IdBNE];
                    g[iC].pf2[IdOUT]=bne.pf2[IdBNE];
                    g[iC].ieqb[IdOUT]=ieqb;
                    g[iC].c1[IdOUT]=mp.cgrid[ic1];
                    g[iC].c2[IdOUT]=mp.cgrid[ic2];
                    g[iC].eqbtype[IdOUT]=(save==1?0:1); //0 = pure, 1 = mixed
                    g[iC].p1[IdOUT]=pstar[j][1]; //because iF was for the opponet few lines above
                    g[iC].p2[IdOUT]=pstar[i][0];
                    g[iC].v11[IdOUT]=br[iC].C0[0][brOUT]+br[iC].C1[0][brOUT]*g[iC].p2[IdOUT];
                    g[iC].v21[IdOUT]=br[iC].C0[1][brOUT]+br[iC].C1[1][brOUT]*g[iC].p1[IdOUT];
                    //simplify
                    A[0]=br[iC].A0[0][brOUT]+br[iC].A1[0][brOUT]*pstar[i][0];
                    B[0]=br[iC].B1[brOUT]*(1-pstar[i][0]);
                    A[1]=br[iC].A0[1][brOUT]+br[iC].A1[1][brOUT]*pstar[j][1];
                    B[1]=br[iC].B1[brOUT]*(1-pstar[j][1]);
                    if (save==1) { //pure strategy
                        //simplify again
                        v0_1[0]=A[0]+B[0]*g[iC].v11[IdOUT];
                        v0_1[1]=A[1]+B[1]*g[iC].v21[IdOUT];
                        for (iF=0;iF<2;iF++) {
                            v0_0[iF]=A[iF]/(1-B[iF]);
                        }
                        g[iC].v10[IdOUT]=v0_1[0]*pbr[0]+v0_0[0]*(1-pbr[0]);
                        g[iC].v20[IdOUT]=v0_1[1]*pbr[1]+v0_0[1]*(1-pbr[1]);
                    }
                    else { //mixed strategy
                        g[iC].v10[IdOUT]=A[0]+B[0]*g[iC].v11[IdOUT];
                        g[iC].v20[IdOUT]=A[1]+B[1]*g[iC].v21[IdOUT];
                    }
                    //reset the save flag
                    save=0;
                }
            }
        }
        //the number of found equilibria
        neqb=ieqb+1;
        g[iC].neqb[id(nC,ic1,ic2,iC,0)]=neqb;
        //mark selected equilibrium
        for (ieqb=0;ieqb<neqb;ieqb++) {
            g[iC].seleqb[id(nC,ic1,ic2,iC,ieqb)]=(ieqb==f_ers(&mp, iC,ic1,ic2,neqb,eqstring)?1.0:0.0);
        }
    }
    else
    {   //eta>0
        //start using brc
        brc.eta=mp.eta;
        brc.B1=br[iC].B1[brOUT];
        for (iF=0;iF<2;iF++) 
        {   //for both functions
            // Coefficients in polynomials short notation (Fedor agreed to this at 15:08 1.12.2011)
            brc.a0[iF]=br[iC].D0[iF][brOUT];
            brc.a1[iF]=br[iC].D1[iF][brOUT];
            brc.a2[iF]=br[iC].D2[iF][brOUT];
        }
        brc.iF=0;//firm 1 (and we will only need iF==1 as seen below :) 
        //iterate over all four root combinations
        //r1 and r2 indicate root combinations
        //and find fixed points
if (ic1==1 && ic2==5) VERBOSE=0;
else VERBOSE=0;
        //initialize ieqb
        ieqb=0;
        for (ir=0;ir<4;ir++)
        {
            brc.r1=r1[ir];
            brc.r2=r2[ir];
            //find the brackets for the segment
            domain[0]=0+mp.ctol;
            domain[1]=1-mp.ctol;
            s_domainbr2(&brc,&mp,domain);
            //if domain is nothing skip this root combination
            if (domain[1]-domain[0]<mp.ctol*2) continue;
if (VERBOSE) printf("iC=%d ic1=%d ic2=%d r1=%d r2=%d domain=[%f,%f]\n",iC,ic1,ic2,(int)r1[ir],(int)r2[ir],domain[0],domain[1]);
            if (r1[ir]*r2[ir]>0)
            {   //this is increasing segment
                //call succ.aprox. algorithm to find possibly multiple
                //equilibria on upward sloping section
                findeqb_sa(&brc,&mp,&mv,domain); //output in mv.eqbinfo
                //if there are any equilibria
                if (mv.neqbinfo>0) {
                    //have to decode the eqbinfo matrix (mv.neqbinfo has number of rows)
if (VERBOSE)
{   //printing
    printf("neqbinfo=%d, eqbinfo:\n",mv.neqbinfo);
    for (i=0; i<5; i++) { // rows
        for (j=0; (j<MAXEQBINFO); j++) { // cols
            printf("%1.6f ", mv.eqbinfo[i+j*5]);
        }
        printf("\n");
    }
    printf("\n");
}
                    //the number of equilibria
                    neqb=2*mv.neqbinfo-3;//normal case
                    //any equilibrium between last and next to last column??
                    if (fabs(mv.eqbinfo[5*(mv.neqbinfo-1)+2]-mv.eqbinfo[5*(mv.neqbinfo-1)+3])<2*mp.ctol) neqb--;
                    //any equilibrium between first and second column??
                    if (fabs(mv.eqbinfo[5*0+2]-mv.eqbinfo[5*0+3])<2*mp.ctol) neqb--;
                    //count equilibria from the last
                    ieqb_offset=neqb-1;
                    //iterate backward over columns in eqbinfo (skip first one)
                    for (i=mv.neqbinfo-1;i>0;i--)
                    {   
                        for (i1=0;i1<2;i1++)
                        {   //possibly two eqb in each column
                            if (i1==0)
                            {   //stable fixed point = stable equilibria
                                if (i==mv.neqbinfo-1) continue;//skip the last column
                                fxpoint=mv.eqbinfo[5*i+0];//first element is fixed point
                            }
                            else
                            {   //unstable fixed point = unstabe equilibria
                                if (neqb==2*mv.neqbinfo-4 && i==mv.neqbinfo-1) continue;//skip last between column in needed
                                if (neqb==2*mv.neqbinfo-5 && i==mv.neqbinfo-1) continue;//skip last between column in needed
                                if (neqb==2*mv.neqbinfo-5 && i==1) continue;//skip first between column in needed
                                fxpoint=(mv.eqbinfo[5*i+2]+mv.eqbinfo[5*(i-1)+3])/2;//first element is fixed point
                            }
if (VERBOSE) printf (">>>>> ieqb=%d fxpoint = %f \n",ieqb+ieqb_offset,fxpoint);
                            //find firm 2 probs
                            

                            //save the equilibrium
                            //ieqb+ieqb_offset is proper index of equilibrium 
                            g[iC].ic1[IdOUT+ieqb_offset]=ic1;
                            g[iC].ic2[IdOUT+ieqb_offset]=ic2;
                            g[iC].pf1[IdOUT+ieqb_offset]=bne.pf1[IdBNE];
                            g[iC].pf2[IdOUT+ieqb_offset]=bne.pf2[IdBNE];
                            g[iC].ieqb[IdOUT+ieqb_offset]=ieqb+ieqb_offset;
                            g[iC].c1[IdOUT+ieqb_offset]=mp.cgrid[ic1];
                            g[iC].c2[IdOUT+ieqb_offset]=mp.cgrid[ic2];
                            g[iC].eqbtype[IdOUT+ieqb_offset]=1; //0 = pure, 1 = mixed
                            g[iC].p1[IdOUT+ieqb_offset]=fxpoint;
                            g[iC].p2[IdOUT+ieqb_offset]=f_invbr(fxpoint,&brc);
                            g[iC].v11[IdOUT+ieqb_offset]=br[iC].C0[0][brOUT]+br[iC].C1[0][brOUT]*g[iC].p2[IdOUT+ieqb_offset];
                            g[iC].v21[IdOUT+ieqb_offset]=br[iC].C0[1][brOUT]+br[iC].C1[1][brOUT]*g[iC].p1[IdOUT+ieqb_offset];
                            g[iC].v10[IdOUT+ieqb_offset]=g[iC].v11[IdOUT+ieqb_offset]+mp.eta*log((1-g[iC].p1[IdOUT+ieqb_offset])/g[iC].p1[IdOUT+ieqb_offset]);
                            g[iC].v20[IdOUT+ieqb_offset]=g[iC].v21[IdOUT+ieqb_offset]+mp.eta*log((1-g[iC].p2[IdOUT+ieqb_offset])/g[iC].p2[IdOUT+ieqb_offset]);
                            //update the offset
                            ieqb_offset--;
                            if (ieqb_offset<0) break; //should never happen!
                        }//i1
                    }
                    //adjust ieqb
                    ieqb+=neqb;
                }
            }
            else 
            {   //this is decreasing segment
                //call bracketing algorithm to find one possible equilibrium
                fxpoint=s_findeqb_bracketing(&brc,&mp,domain);
                if (fxpoint>=0 && fxpoint<=1)
                {
if (VERBOSE) printf("Euilibrium: %f, ieqb=%d\n", fxpoint,ieqb);


                    //save the equilibrium
                    g[iC].ic1[IdOUT]=ic1;
                    g[iC].ic2[IdOUT]=ic2;
                    g[iC].pf1[IdOUT]=bne.pf1[IdBNE];
                    g[iC].pf2[IdOUT]=bne.pf2[IdBNE];
                    g[iC].ieqb[IdOUT]=ieqb;
                    g[iC].c1[IdOUT]=mp.cgrid[ic1];
                    g[iC].c2[IdOUT]=mp.cgrid[ic2];
                    g[iC].eqbtype[IdOUT]=1; //0 = pure, 1 = mixed
                    g[iC].p1[IdOUT]=fxpoint;
                    g[iC].p2[IdOUT]=f_invbr(fxpoint,&brc);
                    g[iC].v11[IdOUT]=br[iC].C0[0][brOUT]+br[iC].C1[0][brOUT]*g[iC].p2[IdOUT];
                    g[iC].v21[IdOUT]=br[iC].C0[1][brOUT]+br[iC].C1[1][brOUT]*g[iC].p1[IdOUT];
                    g[iC].v10[IdOUT]=g[iC].v11[IdOUT]+mp.eta*log((1-g[iC].p1[IdOUT])/g[iC].p1[IdOUT]);
                    g[iC].v20[IdOUT]=g[iC].v21[IdOUT]+mp.eta*log((1-g[iC].p2[IdOUT])/g[iC].p2[IdOUT]);
                    //update ieqb
                    ieqb++;
                }
            }
        }
        //quilibrium selection rule
        //the number of found equilibria
        neqb=ieqb;
        //save neqb to g structure
        g[iC].neqb[id(nC,ic1,ic2,iC,0)]=neqb;
if (neqb==0)
{
printf("iC=%d ic1=%d ic2=%d has no equilibria, decrease ctol!\n",iC,ic1,ic2);
}        
        //mark selected equilibrium
        for (ieqb=0;ieqb<neqb;ieqb++) {
            g[iC].seleqb[id(nC,ic1,ic2,iC,ieqb)]=(ieqb==f_ers(&mp, iC,ic1,ic2,neqb,eqstring)?1.0:0.0);
        }
        
    }
}




//int sa(double* x1, double x0,MPstruct *mp, double l, double r, double lo, double hi) {
int sa(double* x1, double x0, double l, double r,Brcstruct* brc,MPstruct *mp) {
    /*--------------------------------------------------------------------------------------
     * sa() returns domain indicator and fills out x1
     * successive approximations on testfun
     *--------------------------------------------------------------------------------------*/
    int retcode, i;
    double tol, eqb;
/*
    printf("Search for stable equilibrium\n\n");
    
    printf(" i  m          invbr2(m)   tol\n");
*/
    for (i=1; i<=mp[0].maxit; i++){
        //func call
        f_invbr2(x1,x0,brc);

        tol=fabs(x0-x1[0]);
//        printf("%3d %1.8f %1.8f %1.8f\n", i, x0, x1[0], tol);
        if (tol<(mp[0].ctol/10)) {
            eqb=x1[0];
            retcode=0;
//            printf("Equilibrium x1=%g found after %d iterations, tol=%g\n\n", x1[0], i, tol);
            break;
        }
        
        if (x1[0]>=r){
//            printf("above upper bound, r=%g, x1=%g\n\n", r, x1[0]);
            retcode=1;
            break;
        }
        if (x1[0]<=l){
            retcode=-1;
//            printf("below lower bound, l=%g x1=%g\n\n", l, x1[0]);
            break;
        }
        x0=x1[0];
        
    }
    return retcode;
    
}  // end of sa

int comparefct(const void* a, const void* b) {
    double* da = (double*)a;
    double* db = (double*)b;
    
    int diff1 = (da[0] > db[0]) - (da[0] < db[0]);
    if (diff1 != 0) return diff1;
    
    return (da[1] > db[1]) - (da[1] < db[1]);
    
}


//void findeqb(double* eqbinfo) {
void findeqb_sa(Brcstruct* brc,MPstruct *mp,MVstruct *mv,double* domain) {
/*----------------------------------------------------------------------------------------------------------------
 *  This function finds multiple fixed points on the upward sloping sections
    by successive approximations and bracketing
      OUTPUTS:
        mv.eqbinfo:
     rows in eqbinfo: eqb, i, l, r, stable/unstable(1/0);
     cols in eqbinfo: equilibrium info for each equilibrium
        mv.neqbinfo is the number of relevant cols!
     
**----------------------------------------------------------------------------------------------------------------*/
    double lo, hi, l, r, m, x1_m, tol;
    int i, j, k, it, retcode;
    int neqb, ieqb;

    //initial search domain
    lo=domain[0];
    hi=domain[1];
    
// 1: Initialize eqbinfo with NaN and endpoints
    for (j=0;j<(5*MAXEQBINFO);mv[0].eqbinfo[j++]=99); //initialize with >1 numbers
    mv[0].eqbinfo[0+0*5]=lo; //fixed point (or domain borders)
    mv[0].eqbinfo[1+0*5]=1;  //index
    mv[0].eqbinfo[2+0*5]=lo; //lower bound
    mv[0].eqbinfo[3+0*5]=hi; //upper bound
                             //next col is tol
    mv[0].eqbinfo[0+1*5]=hi;
    mv[0].eqbinfo[1+1*5]=2; 
    mv[0].eqbinfo[2+1*5]=lo;
    mv[0].eqbinfo[3+1*5]=hi;
    neqb=2;//initially 2 colums in eqbinfo matrix
    
// 2: start fixed point algorithm
    j=1;
    l=lo;
    r=hi;

    //check for potential no solutions
    //only do the hi point now, lo inside the loop
    if (sa(&x1_m, hi-mp[0].ctol, lo, hi, brc, mp) == -1)
    {   //no eqb, return
        mv[0].neqbinfo=0;
        return;
    }

    //iterate
    for (it=1;it<=mp[0].maxit;it++) {
        //m is next point for function evaluation
        if (it==1){
            m=lo+mp[0].ctol;
        }
        else if (it==2) {
            m=hi-mp[0].ctol;
        }
        else {
            m=(l+r)/2;
        }
        
//        printf("search in interval [l,r]=[%f, %f] starting from m=%f \n\n", l, r, m);
        
        // search for stable equilibira
        retcode=sa(&x1_m, m, l, r, brc, mp);

        //check for potential no solutions
        if (it==1 && retcode==1) 
        {   //started from lo - check for retcode==1
            mv[0].neqbinfo=0;
            return;
        }

        // if equilibirum is found by suucessive approximations
        // Jump to next column to eqbinfo
        if (retcode == 0) {
            neqb=neqb+1;//add a col in eqbinfo matrix
            if (neqb > MAXEQBINFO) mexErrMsgTxt("Error in analytical_sm.c: located equilibria exceed MAXEQBINFO");
            ieqb=neqb-1;//base 0 index
            mv[0].eqbinfo[0 + ieqb*5]=x1_m;
            mv[0].eqbinfo[1 + ieqb*5]=neqb;
            mv[0].eqbinfo[2 + ieqb*5]=l;
            mv[0].eqbinfo[3 + ieqb*5]=r;
            
            //sort by fxpoint and eqb index            
            qsort(mv[0].eqbinfo, MAXEQBINFO, 5*sizeof(double), comparefct);

            //find the index of last eqb after sorting
            for (i=0; i<neqb; i++) { // rows
                if ((int) mv[0].eqbinfo[1+i*5]==neqb) {
                    j=i;
                    break;
                }
            }            
        }
        else {
            x1_m=m;
        }

/*
        //print eqbinfo
        for (i=0; i<5; i++) { // rows
            for (k=0; (k<MAXEQBINFO); k++) { // cols
                if (mv[0].eqbinfo[0+k*5]<99) printf("%8.4f ", mv[0].eqbinfo[i+k*5]);
            }
            printf("\n");
        }
        printf("\n");
*/
        //updating the bounds
        for (k=0;k<neqb;k++) {
            if ((mv[0].eqbinfo[0+k*5]>MAX(x1_m, m)) && (retcode<=0)) {
                mv[0].eqbinfo[2+k*5]=MAX(MAX(x1_m, m), mv[0].eqbinfo[2+k*5]);  // update l if improvement
            }
            else if ((mv[0].eqbinfo[0+k*5]<MIN(x1_m, m)) && (retcode>=0)) {
                mv[0].eqbinfo[3+k*5]=MIN(MIN(x1_m, m), mv[0].eqbinfo[3+k*5]);  // update r if improvement
            }
        }

        //special updates to border columns
        if ((retcode == -1) && (j==0)) { // below l (during left search)(j==1)
            mv[0].eqbinfo[2+j*5]=MAX(m, mv[0].eqbinfo[2+j*5]);
        }
        else if ((retcode ==  1) &&  (j==neqb-1)) {  // above r (during right search)
            mv[0].eqbinfo[3+j*5]=MIN(m, mv[0].eqbinfo[3+j*5]);
        }
        
        //update tolerance
        for (k=1; k<neqb; k++) {
            if (k==0){
            //    mv[0].eqbinfo[4+k*5]=mv[0].eqbinfo[3+k*5]-mv[0].eqbinfo[2+k*5];
            }
            else {
                mv[0].eqbinfo[4+k*5]=mv[0].eqbinfo[3+(k-1)*5]-mv[0].eqbinfo[2+k*5];
            }
        }

        //where next to disect
        tol=0;
        for (k=1; k<neqb; k++) {
            if (tol<mv[0].eqbinfo[4+k*5]) {
                tol=mv[0].eqbinfo[4+k*5];
                j=k; //will disect this interval
                //printf("maxtol=%1.2f", tol);
            }
        }

        //new l and r
        //ctol*5*2 XXXXX
        if (j==0){
            l=mv[0].eqbinfo[2+j*5];
            r=mv[0].eqbinfo[3+j*5]-5*mp[0].ctol*2;
        }
        else if (j==neqb-1){
            l=mv[0].eqbinfo[2+j*5]+5*mp[0].ctol*2;
            r=mv[0].eqbinfo[3+j*5];
        }
        else {
            l=mv[0].eqbinfo[2+j*5]+5*mp[0].ctol*2;
            r=mv[0].eqbinfo[3+(j-1)*5]-5*mp[0].ctol*2;
        }

        //stopping rule
        if (r-l<10*mp[0].ctol){
            break;
        }
    }

    //update the number of equilibria
    mv[0].neqbinfo=neqb;
}


int f_invbr2(double* invbr2, double p_i, Brcstruct* brc) {
/*----------------------------------------------------------------------------------------------------------------
 *  Inverse second order best response function
      OUTPUTS:
      invbr2: inverse second order best response
      returns the indicator of being outside of domain (-1 or +1) or 0 if p_i is inside the domain
**----------------------------------------------------------------------------------------------------------------*/
    double   D, p_j;
    int signinvbr2, jF,iF; 
    double a0,a1,a2;
    //opponent
    iF=brc[0].iF;
    jF=1-brc[0].iF; 
    //rename the coefficients for easy referencing
    a0=brc[0].a0[iF]-brc[0].eta*brc[0].B1*log(p_i)+brc[0].eta*log(p_i/(1-p_i));  
    a1=brc[0].a1[iF]+brc[0].eta*brc[0].B1*log(p_i);
    a2=brc[0].a2[iF];
/*
if (VERBOSE)
{
    printf("iF=%d a0=%f a1=%f a2=%f   jF=%d a0=%f a1=%f a2=%f \n",iF,a0,a1,a2,jF,
        brc[0].a0[jF]-brc[0].eta*brc[0].B1*log(p_i)+brc[0].eta*log(p_i/(1-p_i)),
        brc[0].a1[jF]+brc[0].eta*brc[0].B1*log(p_i),
        brc[0].a2[jF]);
}
*/
    //first discriminant
    D=a1*a1-4*a0*a2;
    signinvbr2=1;//if no roots p_i is too high
    if (D>=0) {
        signinvbr2=-1; //if p_j ouside of unit interval p_i is too low
        p_j=-(1/(2*a2))*(a1+brc[0].r1*sqrt(D));
        if ((p_j)>0 && (p_j)<1) {
            signinvbr2=(int) brc[0].r1; //if no roots, p_i is too low for positive root and too high for negative root
            a0=brc[0].a0[jF]-brc[0].eta*brc[0].B1*log(p_j)+brc[0].eta*log(p_j/(1-p_j));  
            a1=brc[0].a1[jF]+brc[0].eta*brc[0].B1*log(p_j);
            a2=brc[0].a2[jF];
            //second discriminant
            D=a1*a1-4*a0*a2;
            if (D>=0) {
                signinvbr2=-(int)brc[0].r1; //if outside of unit interval p_i is too low for negative root and too high otherwise
                invbr2[0]=-(1/(2*a2))*(a1+brc[0].r2*sqrt(D));
                if ((invbr2[0])>0 && (invbr2[0])<1) { // p_i_min <= p_i =< p_i_max
                    signinvbr2=0; 
                }
                else invbr2[0]=mxGetNaN();
            }
            else invbr2[0]=mxGetNaN(); 
        }
        else invbr2[0]=mxGetNaN();
    }
    else invbr2[0]=mxGetNaN();
    return signinvbr2;
}

double f_invbr(double p_i, Brcstruct* brc) {
/*----------------------------------------------------------------------------------------------------------------
 *  Inverse best response function (first order) for firm iF
    Returns the choice prob of opponent firm 1-iF if own best response is p_i
    (uses coef from firm iF)
**----------------------------------------------------------------------------------------------------------------*/
    double   D;
    double a0,a1,a2;
    int iF; 
    iF=brc[0].iF;
    //rename the coefficients for easy referencing
    a0=brc[0].a0[iF]-brc[0].eta*brc[0].B1*log(p_i)+brc[0].eta*log(p_i/(1-p_i));  
    a1=brc[0].a1[iF]+brc[0].eta*brc[0].B1*log(p_i);
    a2=brc[0].a2[iF];

    //discriminant
    D=a1*a1-4*a0*a2;
    if (D>=0) {
        return -(1/(2*a2))*(a1+brc[0].r1*sqrt(D));
    }
    else return mxGetNaN();
}

//void s_findminmaxbr2(double pimin,double pimax, double invbr2_min,double invbr2_max, 
//                     double a0[2], double a1[2], double a2[2], double B1, 
//                     MPstruct *mp, int iF, int ir) {
void s_domainbr2(Brcstruct *brc,MPstruct *mp,double* domain) {
/*----------------------------------------------------------------------------------------------------------------
 *  Min and max of the support of the inverse second order best response function (by segment)
 *  i.e. domain of the second order best response function (by segment)
    INPUT: search region = domain
    OUTPUT: domain
**----------------------------------------------------------------------------------------------------------------*/
        double u, l, tmp;
        int it;

        l=domain[0];
        u=domain[1];
        //search for max
        if (f_invbr2(&tmp,u,brc)==1)
        {   //search only if there is space above, otherwise keep the upper bound
             it=0;
            while (fabs(u-l)>mp[0].ctol && it<=mp[0].maxit)
            {
                if (f_invbr2(&tmp,(u+l)/2,brc)>0) u=(u+l)/2;
                else l=(u+l)/2;
                it++;
            }
            domain[1]=l;
        }

        l=domain[0];
        u=domain[1];
        //search for min
        if (f_invbr2(&tmp,l,brc)==-1)
        {   //search only if there is space below, otherwise keep the lower bound
            it=0;
            while (fabs(u-l)>mp[0].ctol && it<=mp[0].maxit)
            {
                if (f_invbr2(&tmp,(u+l)/2,brc)>-1) u=(u+l)/2;
                else l=(u+l)/2;
                it++;
            }
            domain[0]=u;
        }
}

double s_findeqb_bracketing(Brcstruct *brc,MPstruct *mp,double* domain) {
    /*This function compute equilibrium by searching for fixed point on the
     * inverse second order best reponse function.
     * Only called for downward sloping segments!
     */
        double u, l, tmp;
        int it;

        l=domain[0];
        u=domain[1];
        //check if there is fixed point at all at this segment
        //check the lower bound
        f_invbr2(&tmp,l,brc);        
        if (tmp-l<0) return -1; //no equilibria here
        //check the upper bound
        f_invbr2(&tmp,u,brc);
        if (tmp-u>0) return -1; //no equilibria here

        //search for fixed point
        it=0;
        while (fabs(u-l)>mp[0].ctol && it<=mp[0].maxit)
        {
//printf(">>> it=%d u=%f l=%f invbr2(mid)=%d %f\n",it,u,l,f_invbr2(&tmp,(u+l)/2,brc),tmp);
            if (f_invbr2(&tmp,(u+l)/2,brc),tmp-((u+l)/2)<0) u=(u+l)/2;
            else l=(u+l)/2;
            it++;
        }
        return (u+l)/2;
//printf(">>> it=%d\n",it-1);

}


