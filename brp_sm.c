//         Fedor Iskhakov, University Technology Sidney
//         John Rust, University of Maryland
//         Bertel Schjerning, University of Copenhagen
//         March 2012

// Method for looking for equilibria: EQBMETHOD_GRID, EQBMETHOD_VF
#define EQBMETHOD_GRID

// Switch for drawing equilibrium graph EQB_VISUAL
// #define EQB_VISUAL

double f_brp1_sm(double p2, MPstruct* mp, Brc_amstruct* brc) {
    brc[0].v1[0]=brc[0].pf[0]-brc[0].kv+mp[0].df*(1-brc[0].pc)*brc[0].logsumK     +p2*brc[0].x0[0]+(1-p2)*brc[0].x1[0];
    brc[0].v0[0]=f_solvevf(brc[0].A0[0]*p2+brc[0].A1[0]*(1-p2), brc[0].B[0]*(1-p2), brc[0].v1[0], mp[0].eta,   mp);
    brc[0].p[0]=f_logit(2, mp[0].eta, NULL, brc[0].v1[0], brc[0].v0[0]);
    return brc[0].p[0];
}

double f_brp2_sm(double p1, MPstruct* mp, Brc_amstruct* brc) {
    brc[0].v1[1]=brc[0].pf[1]-brc[0].kv+mp[0].df*(1-brc[0].pc)*brc[0].logsumK     +p1*brc[0].x0[1]+(1-p1)*brc[0].x1[1];
    brc[0].v0[1]=f_solvevf(brc[0].A0[1]*p1+brc[0].A0[1]*(1-p1), brc[0].B[1]*(1-p1), brc[0].v1[1], mp[0].eta,   mp);
    brc[0].p[1]=f_logit(2, mp[0].eta, NULL, brc[0].v1[1], brc[0].v0[1]);
    return brc[0].p[1];
}

double  f_br2p1_sm(double p1, int ieqb, int eqbtype, MPstruct* mp, Brc_amstruct* brc, Gamestruct* gc, int ic1, int ic2, int iC, int nC) {
    double p2;
    // FIRM 2's response to the perception that firm 1 will invest with probability p1
    p2=f_brp2_sm(p1, mp, brc);     

    //  FIRM 1's second order response to the perception that it will invest with probability p1
    brc[0].p[0]=f_brp1_sm(p2, mp, brc);
    
    // save results in g structure whenever ieqb>0 
    if (ieqb>=0){
        gc[0].ieqb[IdOUT]=ieqb;
        gc[0].neqb[id(nC,ic1,ic2,iC,0)]++;
        gc[0].eqbtype[IdOUT]=eqbtype;
        gc[0].v10[IdOUT]=brc[0].v0[0];
        gc[0].v11[IdOUT]=brc[0].v1[0];
        gc[0].v20[IdOUT]=brc[0].v0[1];
        gc[0].v21[IdOUT]=brc[0].v1[1];
        gc[0].p1[IdOUT]=brc[0].p[0];
        gc[0].p2[IdOUT]=p2;
    }
    return (brc[0].p[0]);
}

double  f_br1p2_sm(double p2, int ieqb, int eqbtype, MPstruct* mp, Brc_amstruct* brc, Gamestruct* gc, int ic1, int ic2, int iC, int nC) {
    double p1;
    // FIRM 2's response to the perception that firm 1 will invest with probability p1
    p1=f_brp1_sm(p2, mp, brc);     

    //  FIRM 1's second order response to the perception that it will invest with probability p1
    brc[0].p[1]=f_brp2_sm(p1, mp, brc);
    
    // save results in g structure whenever ieqb>0 
    if (ieqb>=0){
        gc[0].ieqb[IdOUT]=ieqb;
        gc[0].neqb[id(nC,ic1,ic2,iC,0)]++;
        gc[0].eqbtype[IdOUT]=eqbtype;
        gc[0].v10[IdOUT]=brc[0].v0[0];
        gc[0].v11[IdOUT]=brc[0].v1[0];
        gc[0].v20[IdOUT]=brc[0].v0[1];
        gc[0].v21[IdOUT]=brc[0].v1[1];
        gc[0].p1[IdOUT]=p1;
        gc[0].p2[IdOUT]=brc[0].p[1];
    }
    return (brc[0].p[1]);
}

