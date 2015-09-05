//         Fedor Iskhakov, University Technology Sydney
//         John Rust, University of Maryland
//         Bertel Schjerning, University of Copenhagen
//         March 2012

// Method for looking for equilibria: EQBMETHOD_GRID, EQBMETHOD_VF
#define EQBMETHOD_GRID

// All possible equilibrium selection rules
// MUST BE THE SAME AS IN SETUP.M
#define ESR_MixedStr             1
#define ESR_FirstInvests         2
#define ESR_HighCostInvestsMix   3
#define ESR_HighCostInvests      4
#define ESR1                     5
#define ESRAntiLeapfrogging      6
#define ESRLeapfrogging          7
#define ESRstring               99  

int iK (int iC, int ic1, int ic2,int nC) {
    /*----------------------------------------------------------------------------------------------------------------
     *          Index in the eqstring
     **----------------------------------------------------------------------------------------------------------------*/
    //Calculates the index in the eqstring that corresponds to the given 
    //iC,ic1,ic2 AND preserves the order of dependence of solutions so
    //that lexicographical order in eqstrings allows not to recalculate
    //the solution of same and lower layers (which lie in the higher digits
    //of eqstring)
    int K=0;
    int iCp;
    
    //skip blocks by iC
    iCp=nC-iC-1;
    //sum of squares from 1 to iCp
    //while (iCp>0) {K+=iCp*iCp;iCp--;}
    K+=iCp*(iCp+1)*(2*iCp+1)/6;
    
    //interior
    if (ic1!=iC && ic2!=iC) return K+=(nC-iC-1)*(ic1-iC-1) + (ic2-iC-1);
    K+=(nC-iC-1)*(nC-iC-1);
    
    //edges    
    if (ic1!=iC && ic2==iC) return K+=ic1-iC-1;
    K+=(nC-iC-1);
    if (ic1==iC && ic2!=iC) return K+=ic2-iC-1;
    K+=(nC-iC-1);
    
    //corner - one element, so nothing to add   
    return K;
}    

void iKinv (int K,int nC,int* iC,int* ic1,int* ic2) {
    /*----------------------------------------------------------------------------------------------------------------
     *          ic1,ic2,iC for a given index in eqsting
     **----------------------------------------------------------------------------------------------------------------*/
    //Calculates the triple (ic1,ic2,iC) from the index in the eqstrin
    
    iC[0]=nC-1;
    //skip blocks by iC
    while (K>=(nC-iC[0])*(nC-iC[0]))
    {   K-=(nC-iC[0])*(nC-iC[0]);
        iC[0]--;
    }
    //iC has right value on exit from while

    //interior
    if (K<(nC-iC[0]-1)*(nC-iC[0]-1))
    {
        ic1[0]=K/(nC-iC[0]-1) +iC[0]+1;
        ic2[0]=K%(nC-iC[0]-1) +iC[0]+1;
        return;
    }
    K-=(nC-iC[0]-1)*(nC-iC[0]-1);

    //edges
    if (K<nC-iC[0]-1)
    {
        ic2[0]=iC[0];
        ic1[0]=K +iC[0]+1;
        return;
    }
    K-=(nC-iC[0]-1);
    if (K<nC-iC[0]-1)
    {
        ic2[0]=K +iC[0]+1;
        ic1[0]=iC[0];
        return;
    }
    
    //the corner
    ic2[0]=iC[0];
    ic1[0]=iC[0];
    return;
}    

//just for testing XXX
void prepareEqstring(int *eqstring, int nC) {
    //create the ewqstring using current rule ESR from modelparts.c
    //assumes eqstring has enough space!!!!
    //CAN ONLY BE USED IN GAMES ALREDY STUDIED 
    int ic1,ic2,iC;
    int k,ic1n,ic2n,iCn;
    for (iC=0;iC<nC;iC++)
    {
        for (ic1=iC;ic1<nC;ic1++)
        {
            for (ic2=iC;ic2<nC;ic2++)
            {
                //first equilibrium
                //if (ic1==iC || ic2==iC) eqstring[iK(iC,ic1,ic2,nC)]=0;
                //else eqstring[iK(iC,ic1,ic2,nC)]=1;
                k=iK(iC,ic1,ic2,nC);
//                printf("iC=%d ic1=%d ic2=%d --> K=%d \n",iC,ic1,ic2,k);
                iKinv(k,nC,&iCn,&ic1n,&ic2n);
                if (k!=iK(iCn,ic1n,ic2n,nC)) printf("FAIL!!!!");
//                printf("iC=%d ic1=%d ic2=%d <-- K=%d \n\n",iCn,ic1n,ic2n,k);
            }
        }
    }
}

int f_ers(MPstruct *mp, int iC, int ic1, int ic2, int neqb, int *eqstr) {
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
//
    int res;
    switch (mp[0].esr) {
        case ESRstring:
            return eqstr[iK(iC,ic1,ic2,mp[0].nC)];
            //don't do anything else for ESRstring for speed    
        case ESR1:
            res=0;
            break;    
        case ESR_HighCostInvests:   
            if ((ic1 == iC) || (ic2 == iC)) res=0;
            else if (neqb == 1) res=0;
            // The firm with the highest cost gets to invest and in case of a tie, firm 1 gets to invest
            else if (ic1 < ic2) res=0;          // firm 2 gets to invest
            else if (ic1 > ic2) res=neqb-1;     // firm 1 gets to invest
            else if (ic1 == ic2) res=neqb-1;    // firm 1 gets to invest in case of tie
            break;    
        case ESR_HighCostInvestsMix:   
//             if ((ic1 == iC) || (ic2 == iC)) res=0;
            if (neqb == 1) {
                res=0;
            }
            else {
                // The firm with the highest cost gets to invest and in case of a tie, firm 1 gets to invest
                if (ic1 < ic2) res=0;          // firm 2 gets to invest
                if (ic1 > ic2) res=neqb-1;     // firm 1 gets to invest
                if (ic1 == ic2) {
                    if (neqb==3) res=1;
                    if (neqb==5) res=2;
                }
            }
            break;    
        case ESR_MixedStr:   
            if ((ic1 == iC) || (ic2 == iC)) res=0;
            else if (neqb == 1) res=0;
            else if (neqb==3) res=1;
            else if (neqb==5) res=2;
            else res=0;
            break;    
        case ESR_FirstInvests:   
            if ((ic1 == iC) || (ic2 == iC)) res=0;
            else if (neqb == 1) res=0;
            else res=neqb-1;
            break;    
        case ESRAntiLeapfrogging:    
            if ((ic1 == iC) || (ic2 == iC)) res=0;
            else if (neqb == 1) res=0;
            else if (ic1==ic2) res=1;
            else if (ic1>ic2) res=0;
            else if (ic1<ic2) res=2;
            break;
        case ESRLeapfrogging:    
            if ((ic1 == iC) || (ic2 == iC)) res=0;
            else if (neqb == 1) res=0;
            else if (ic1==ic2) res=1;
            else if (ic1>ic2) res=0;  // firm2 invests
            else if (ic1<ic2) res=2;  // firm1 invests
            break;
        default:
            mexErrMsgTxt("Unknown equilibrium selection rules or no equilibrium selected for given values of state variables!");
    }

    //print out the result:
    //printf("iC=%d ic1=%d ic2=%d neqb=%d K=%d eqb=%d\n",iC,ic1,ic2,neqb,iK(iC,ic1,ic2,mp[0].nC),res);
    //printf("%d %d\n",iK(iC,ic1,ic2,mp[0].nC),res);
    if (mp[0].out4) mp[0].out4[(mp[0].nC*(mp[0].nC+1)*(2*mp[0].nC+1)/6)-1-iK(iC,ic1,ic2,mp[0].nC)]=res;  //if pointer is not NULL record the ESR rule as string

    //return the result
    return res;

}
