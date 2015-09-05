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

global par mp;
%% mp: Structure that hold model parameters

mp.df=0.95;
% per period discount factor for the firms' discounted profits

mp.k1=0.3;
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

mp.eta=0.0;
% scaling parameter for unobserved extreme value investment cost shocks affecting the duopolists'
% decisions on whether to invest in an improved state of the art production facility
% or not. When eta is zero, there are no unobserved investment cost shocks, and each
% firm invests in the state of the art if the expected future discounted profits net of
% the cost of investing is greater than the discounted future profits of continuing to
% produce using the firm's existing production equipment

mp.c_og=5.0;
% the "cost" of the outside good, if one is present in the model

mp.c_tr=0.1;
% coefficient on the current cost of production (representing state of
% current technology) in the probability of a further advancement in
% technology this period, represented by a further lowering of costs,
% as given by a draw from a beta distribution on the interval [c,cmin].


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

par.nC=50;
% number of grid points for the cost variables (same number of grids used for
% each firm)

par.nP=60;
% number of points on the [0,1] interval used for initial grid search to find fixed points in
% seg.m and regs.m.  Used to create p1vec vector below


par.ctol=1e-9;
% convergence tolerance for Newton algorithm

%% OLD and INACTIVE STUFF BELOW THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  OLD STUFF     OLD STUFF      OLD STUFF      OLD STUFF            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
   esr=3;                % equilibrium selection rule:
                         % 0: select lowest probability of (firm 1) investing whenever there are multiple equilibria
                         % 1: select middle equilibrium (mixed strategy for firm 1) whenever there are multiple equilibria
                         % 2: select highest probability of firm 1 investing whenever there are multiple equilibria
                         % 3: select the deterministic leap-frogging equilibrium: 
                         %    when c1 < c2,  firm 1 invests with probability 0, firm 2 invests when c is sufficiently below c1
                         %    when c1 > c2,  firm 2 invests with probability 0, firm 1 invests when c is sufficiently below c2
                         %    when c1 = c2,  firms 1 and 2 play the mixed strategy (middle) equilibrium


p1sp =1.0;               % grid spacing parameter for the p1vec [XXX currently not used]

   c1=5.0;               % marginal cost of production for firm 1 under its initial production technology

   c2=5.0;               % marginal cost of production for firm 2 under its initial production technology


   clb=0.8;              % fraction of the current production cost c representing the largest
                         % possible decrease in costs in any given period


   csp=1.0;              % spacing parameter for the cost grid point: csp=1 results in uniformly
                         % spaces points between [cmin,cmax] whereas csp > 1 results in points that
                         % are more closely spaced near cmin, and csp < 1 results in points that are
                         % more closely spaces near cmax

                         

   beta_a=0.0;           % a parameter of the beta distribution over the interval [c,cmin] where
                         % c is the current best technology marginal cost of production and cmin
                         % is a lower bound on this marginal cost beyond which technological improvements
                         % can no further reduce costs, even in the limit as time goes to infinity.
                         % So the expected cost, given an innovation occurs, will be
                         % cmin+(c-cmin)*beta_a/(beta_a+beta_b)

   beta_b=0.0;           % a parameter of the beta distribution over the interval [c,cmin] where
                         % c is the current best technology marginal cost of production and cmin
                         % is a lower bound on this marginal cost beyond which technological improvements
                         % can no further reduce costs, even in the limit as time goes to infinity.
                         % So the expected cost, given an innovation occurs, will be
                         % cmin+(c-cmin)*beta_a/(beta_a+beta_b)



   nqp=10;               % number of quadrature points used in the numerical integration

   satol=1e-6;           % convergence tolerance for successive approximations algorithm for (v10,v11,v20,v21)


   init1=1;              % sets the order of the backward induction calculations in bellman.c:
                         % init1=1, do the calculations for firm 1 first, then do firm 2
            			 % init1=0, do the calculations for firm 2 first, then do firm 1

   intmeth=1;            % 1 for 3 dimensional multi-linear interpolation, 0 for simplicial interpolation
                         % of the value functions in the successive approximation iterations in dds_c.c

   pr=1;                 % 1 for detailed output on the iterations (useful in debugging)
                         % 0 for only summary output

   debg=0;               % 1 for additional output for use in debugging program
                         % 0 for no debugging output

%}


