#define OTHER_INPUTS , MPstruct* mp, Brc_amstruct* brc, Gamestruct* gc, int ic1, int ic2, int iC, int nC
#define OTHER_INPUTS_CLEAN , mp, brc, gc, ic1, ic2, iC, nC
void s_ffxp(double (*fun)(double p1, int eqnum, int eqtype OTHER_INPUTS) OTHER_INPUTS) {
    /*----------------------------------------------------------------------------------------------------------------
Routine to look for equilibria by solution of the p1=f(p1,...) function on [0,1] interval.
During iterations f is called with command=-1 and should not save any results, just return f(p1,...) as fast as possible.
When equilibrium K=0,1,2,.. is found f is called as f(p*,K,...) where p* is p1 in equilibrium.
Assumptions:
Edit for the f other inputs
void s_ffxp (double (*f)(double p1, int command, <other imputs>), <same other imputs>) {
     *----------------------------------------------------------------------------------------------------------------*/
    int i; //indexes
    const double eps=0.000001; //small number for calculation of "derivatives"
    double x,f,xprev,fprev;
    double brck1, brck2, brck3, brf1, brf2, brf3, tolmeasure, tolmeasureprev, looptest; //variables for bracketing
    int brckmethod, it; //indicator of bracketing method to be used
    int sig,sigprev; //indicator of above =1 or below =0 45 degree line
    int eqbcount=0, eqbtype; //equibilrium count, type (see seg.c or below for description)
    int savelast=0; //indicator for saving the pure strategy pr=1 equilibrium last
#ifdef EQB_VISUAL
mxArray *Mvars[4];
const int Mpoints=1000;
double *Mvar;
char Mstring[500];
/*Plot best responce function*/
for (i=0;i<4-1;i++) {
    Mvars[i]=mxCreateDoubleMatrix(Mpoints,1,mxREAL);
    Mvar=mxGetPr(Mvars[i]);
    for (j=0;j<Mpoints;Mvar[j++]=(double)j/Mpoints);
}
Mvars[3]=mxCreateDoubleMatrix(Mpoints,1,mxREAL);
Mvar=mxGetPr(Mvars[3]);
for (j=0;j<Mpoints;Mvar[j++]=(*fun)((double)j/Mpoints,-1,-1 OTHER_INPUTS_CLEAN));
mexCallMATLAB(0,NULL,4,Mvars,"plot");
sprintf(Mstring,"title('ic1=%d ic2=all ic=%d')",ic1,iC);
mexEvalString(Mstring);
mexEvalString("hold on");
if (ic2==nC-1) {
    mexEvalString("drawnow");
    mexEvalString("input ('Press Enter to continue..');");
    mexEvalString("close gcf;");
}
#endif
#ifdef EQBMETHOD_GRID
// Grid search + bracketing
for (i=0;i<=mp[0].nP;i++) {
    x=(double)i/(double)mp[0].nP;
    f=(*fun)(x,-1,-1 OTHER_INPUTS_CLEAN);
    if (fabs(x-f)<1e-10) {
        //found eqb
        if (eqbcount>=MAXEQB) mexErrMsgTxt("Found more equilibria than allowed by MAXEQB!"); //quit to Matlab if >MAXEQB
        if (fabs(x)<1e-10) { //pure strategy with pr=0
            x=eps;
            f=(*fun)(x,-1,-1 OTHER_INPUTS_CLEAN);
            brckmethod=1;
            if (f<x) {
                eqbtype=1; //pure stable
                sig=0; //as if sigprev was 1
            }
            else {
                eqbtype=0; //pure unstable
                sig=1; //as if sigprev was 0
            }
            (*fun)(0.0,eqbcount,eqbtype OTHER_INPUTS_CLEAN);
            eqbcount++;
            xprev=x;
            fprev=f;
            sigprev=sig;
            continue;
        }
        else if (x>0.0 && x<1.0) { //mixed strategy, accidentally hit the eqb point
            //let the section on mixed eqb handle this
            xprev=x-eps;
            fprev=(*fun)(xprev,-1,-1 OTHER_INPUTS_CLEAN);
            sigprev=(int)(fprev>xprev);
            x=x+eps;
            f=(*fun)(x,-1,-1 OTHER_INPUTS_CLEAN);
            sig=(int)(f>x);
        }
        else if (fabs(x-1.0)<1e-10) { //pure strategy with pr=1
            x=1-eps;
            f=(*fun)(x,-1,-1 OTHER_INPUTS_CLEAN);
            brckmethod=1;
            //save pure strategy equilibrium
            if (f>x) {
                eqbtype=1; //pure stable
                sig=1; //as if sig was sigprev=1
            }
            else {
                eqbtype=0; //pure unstable
                sig=0; //as if sig was sigprev=0
            }
            savelast=eqbtype;//delay saving this equlibrium
            //let mixed strategy check for if sig changed
        }
    }
    else sig=(int)(f>x);
    if (x>0 && sig!=sigprev) {
        //found eqb! (mixed strategy only in this part)
        if (eqbcount>=MAXEQB-(savelast>0?1:0)) mexErrMsgTxt("Found more equilibria than allowed by MAXEQB!"); //quit to Matlab if >MAXEQB
        if (sigprev==1) {eqbtype=3;} else {eqbtype=2;} //eqlb type
        //bracketing
        brck1=xprev;
        brf1=fprev-xprev;
        brck2=x;
        brf2=f-x;
        brck3=(brck1+brck2)/2;
        tolmeasure=fabs(brck2-brck1);
        tolmeasureprev=1.0; //initial values for the "loop" test
        it=0;
        while (tolmeasure>mp[0].ctol && it<mp[0].maxit) {
            it++;
            if (brckmethod==0) brck3=(brf1*brck2-brf2*brck1)/(brf1-brf2); //linear interpolation of the function in question (like successive average)
            else brck3=(brck1+brck2)/2; //binary section (method 1)
            brf3=(*fun)(brck3,-1,-1 OTHER_INPUTS_CLEAN)-brck3;
//printf ("f1(%0.3f)=%0.3f f2(%0.3f)=%0.3f -> f3(%0.3f)=%0.3e (%s)",brck1,brf1,brck2,brf2,brck3,brf3,(brckmethod==0?"linear appx bracketing":"binary sections"));
            if ((brf1>0 && brf3>0) || (brf1<0 && brf3<0)) {
                tolmeasure=fabs(brck3-brck1);
                looptest=(brf3-brf1)*brf2*(brck3-brck2)/((brf2-brf3)*(brf2-brf3)*(brck3-brck1)) + brf2/(brf2-brf3); //see loop test below
                brck1=brck3;
                brf1=brf3;
//printf ("<tol=%0.4f> 1->3\n",tolmeasure);
            }
            else if ((brf2>0 && brf3>0) || (brf2<0 && brf3<0)) {
                tolmeasure=fabs(brck3-brck2);
                looptest=(brf3-brf2)*brf1*(brck3-brck1)/((brf1-brf3)*(brf1-brf3)*(brck3-brck2)) + brf1/(brf1-brf3); //see loop test below
                brck2=brck3;
                brf2=brf3;
//printf ("<tol=%0.4f> 2->3\n",tolmeasure);
            }
            else if (fabs(brf3)<1e-10) {
                tolmeasure=0.0;
//printf ("<tol=%0.4f> stop\n",tolmeasure);
            }
            /* "Loop" test: test for the relative size of the step
             * It can be shown that the limit of the retalive step is approaching
             * f'(x1)*(y0*(x1-x0)/(y0-y1)^2) + y0/(y0-y1) where the x1 is the
             * changing side of the bracket and y1=f(x1), y2=f(x2)
             * Therefore switch to bisections when this limit is anywhere close
             */
            if (fabs(tolmeasure/tolmeasureprev-looptest)<0.5*fabs(looptest)) brckmethod=1;
            tolmeasureprev=tolmeasure;
        }
        if (it>=mp[0].maxit) mexErrMsgTxt("Maximum number of iterations exceeded in s_ffxp!");
        (*fun)(brck3,eqbcount,eqbtype OTHER_INPUTS_CLEAN);
        eqbcount++;
    }
    xprev=x;
    fprev=f;
    sigprev=sig;
    brckmethod=0; //to start next interval clean
}
//save pr=1 pure strategy equilibrium last (for nice sorting w.r.t. pr1)
if (savelast>0) (*fun)(1.0,eqbcount,savelast OTHER_INPUTS_CLEAN);
#endif
#ifdef EQBMETHOD_VF
mexErrMsgTxt("Equilibrium search with value function type iterations is yet to be written!");
/* Ideas for optimization of this routine:
 * Available tools - value function iterations, newton method, bracketing, grid search
 * if reaction function can be inverced unstable equilibria can be found
 * with value function type iterations too
 * Assume safely: odd number of equilibria, first and last eqb is stable
 *
 */
#endif
}


