%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                   %%%%
%%%%    setup.m: setup global parameter structures for the solution    %%%%
%%%%             to the dynamic duopoly problem                        %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                   %%%%
%%%% BY:     Fedor Iskhakov, University Technology Sidney              %%%%
%%%%         John Rust, University of Maryland                         %%%%
%%%%         Bertel Schjerning, University of Copenhagen               %%%%
%%%%                                                                   %%%%
%%%% THIS VERSION: March 2011                                          %%%%
%%%%                                                                   %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global par mp sw;
%% mp: Structure that hold model parameters
mp.dt=1;                 
% length of time intervals. Let this go to zero to compute the equilibria of the limiting
% continuous time game

mp.df=exp(-0.05*mp.dt);
% per period discount factor for the firms' discounted profits

% mp.k1=28.3; %gives multiple eqbria
%mp.k1=18.3;
mp.k1=8.3;

% parameter in the numerator of the function k(c), representing the cost of
% acquiring plant capable of producing at the current state of the art marginal
% cost of production, c

mp.k2=1.0;
% parameter the denominator, the coefficient of c, in the function k(c)

mp.sigma=0.00000;
% extreme value scaling parameter for the effect of unobserved demand shocks affecting
% consumers' choices between the two goods produced the duopolists and the outside good,
% if present. As sigma goes to zero, each consumer simply chooses the good that is
% the cheapest, or the outside good if the "cost" of the outside good is less than the
% Bertrand duopoly prices charged by the two firms

mp.eta=0;
% scaling parameter for unobserved extreme value investment cost shocks affecting the duopolists'
% decisions on whether to invest in an improved state of the art production facility
% or not. When eta is zero, there are no unobserved investment cost shocks, and each
% firm invests in the state of the art if the expected future discounted profits net of
% the cost of investing is greater than the discounted future profits of continuing to
% produce using the firm's existing production equipment

mp.c_og=5.0;
% the "cost" of the outside good, if one is present in the model

mp.c_tr=0.05;
mp.c_tr=-1;
% coefficient on the current cost of production (representing state of
% current technology) in the probability of a further advancement in
% technology this period, represented by a further lowering of costs,
% as given by a draw from a beta distribution on the interval [c,cmin].
% if mp.c_tr<0 technonoly improves deterministically each period until c=0

mp.beta_a=1.8; 
% a parameter of the beta distribution over the interval [c,cmin] where
% c is the current best technology marginal cost of production and cmin
% is a lower bound on this marginal cost beyond which technological improvements
% can no further reduce costs, even in the limit as time goes to infinity.
% So the expected cost, given an innovation occurs, will be
% cmin+(c-cmin)*beta_a/(beta_a+beta_b)

mp.beta_b=0.4;           
% b parameter of the beta distribution over the interval [c,cmin] where
% c is the current best technology marginal cost of production and cmin
% is a lower bound on this marginal cost beyond which technological improvements
% can no further reduce costs, even in the limit as time goes to infinity.
% So the expected cost, given an innovation occurs, will be
% cmin+(c-cmin)*beta_a/(beta_a+beta_b)
             
mp.onestep=1;  
% set to one to replicate results for original model where technological
% improvements can only move to the next grid point down

mp.tpm=[0.5 .5;
        0.5 .5]; 
    
% mp.tpm=[0.2 0.8; 0.6  0.4];  % mover transition probability matrix

% mp.tpm=[0.0 1.0; 1.0  0.0]; % mover transition probability matrix

%% par: Structure that hold parameters that index spaces, algorithms, etc
par.og=0;
% enter 0 for solution to problem with no outside good,
% otherwise enter 1 to include an outside good
% treated as integer in c

par.cmin=0.0;
% lower bound on the state space, the lowest possible value of the
% state of teh art marginal cost of production

par.cmax=5.0;
% upper bound on the state space, the highest possible value of the
% state of teh art marginal cost of production

par.maxit=200;
% maximum number of successive approximations iterations allowed

par.nC=5;
% number of grid points for the cost variables (same number of grids used for
% each firm)

par.nP=100;
% number of points on the [0,1] interval used for initial grid search to find fixed points in
% seg.m and regs.m.  Used to create p1vec vector below

par.ctol=1e-12;
% convergence tolerance for Newton algorithm

par.pti=f_pti(mp, par);
% matrix to hold the technological improvement probabilities 
% at technology state c=cgrid(j), this is a probability distribution p(0),...,p(cgrid(j-1))
% (conditional on tech. improvement)


%% sw: Structure that hold parameters that index equilibrium selection rules
sw.alternate=false; %alternate move or simultanious move game
%%Equilibrium selection rules
% ESR_MixedStr             =1;
% ESR_FirstInvests         =2;
% ESR_HighCostInvestsMix   =3;
% ESR_HighCostInvests      =4;
% ESR1                     =5;
% ESRstring               =99;  

sw.analytical=false; 
% solve the model using analytical solutions when eta>0
% NOTE: Analytical solution only implemented active for the simultanious move
%       game when eta>0
% When eta=0, model is solved always solved using analytical solution


sw.esr=3; 
%equilibrium selection rule to be used: see setup.m or esr.c

sw.esrmax=100000; 
%run N feasible eqstrings (set equal to 192736405 for n=5)

sw.esrstart=0; 
%index of the first eqstring


%names of the colums of output data
gcol={'1: ic1' '2: ic2' '3: ieqb' '4: c1' '5: c2' '6: eqbtype (NA)' '7: v10 firm1 not' '8: v11 firm1 invest' '9: v20 firm2 not' '10: v21 firm2 invest' ...
      '11: p1 firm1 invest' '12: p2 firm2 invest' '13: seleqb' '14: pf1 firm1' '15: pf2 firm2' '16: x10' '17: x11' '18: x20' '19: x21' '20: neqb' ...
      };
gcol_ec={'1: ecv10' '2: ecv11' '3: ecv20' '4: ecv21' '5: ecx10' '6: ecx11' '7: ecx20' '8: ecx21','9: ecm1 total','10: ecm2 total'};


