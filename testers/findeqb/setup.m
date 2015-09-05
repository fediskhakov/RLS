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
