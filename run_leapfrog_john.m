%         Fedor Iskhakov, University Technology Sydney
%         John Rust, Georgetown University
%         Bertel Schjerning, University of Copenhagen
%         Sydney, August 2012


clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPILE the C code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compiler flags:
% MAXEQB - maximum number of equilibria in a stage game, used for space allocation
% PRINTeqbloop - amount of printing in the loop over the lexicographical equilibrium
%          selection rules (ESR). 0=nothing, 1=results, 2=ESR strings
% PRINTeqbstr - whether to print ESR strings for found equilibria (0 or 1)
% more compiler flags can be found in the top of leapfrog.c file
fprintf('Compiling.. ');tic;
mex -largeArrayDims leapfrog.c -DMAXEQB=3 -DPRINTeqbloop=0 -DPRINTeqbstr=0
tc=toc; fprintf('done in %1.10fsec\n',tc);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <<1>> load default parameters
setup

% <<2>> adjust parameters one by one
% The main parameters are:
% > time differential
mp.dt=1;
% > investment cost
mp.k1=10;
mp.k2=1.0;
% > eta
mp.eta=0;
% > technological progress
mp.c_tr=1;
mp.onestep=1;  
mp.beta_a=1.8; 
mp.beta_b=0.4;           
% > transition probability for turns (standard indexing: from in rows, to in colums)
mp.tpm=[0 1;
        1 0]; 
% > number of grid points
par.nC=4;

% <<3>> recalculate depend parameters
%     !!!! dependence on dt needs to be included here !!!!
par.pti=f_pti(mp, par);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOLVE the model ONCE with simple ESR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <<1>> set the model type and ESR
sw.alternate=false;% alternate move or simultanious move game
sw.analytical=true;% use analytical of numerical solution, this settings applies
                   % ONLY for simultanous move, eta>0  
% <<2>> choose ESR
% ESR_MixedStr             =1;
% ESR_FirstInvests         =2;
% ESR_HighCostInvestsMix   =3;
% ESR_HighCostInvests      =4;
% ESR1                     =5;
% ESRstring               =99;
% Precise definitions are in esr.c
sw.esr=5; %equilibrium selection rule to be used: see setup.m or esr.c

% If sw.esr is not ESRstring (sw.esr<99) the model is solve only once, and using
% one of the ESR among those we used before.  ESR is the same for all points of 
% the state space

% <<3>> run the solver
fprintf('Running the solver.. ');tic;
[bne, br, g, eqbstr]=leapfrog(par,mp,sw);
ts=toc;fprintf('done in %1.10fsec\n',ts);

% Output is:
%   per-period payoffs
%   coefficients from analytical solution (ONLY in simultanous move game)
%   game structure g
%   table with equilibrium statistics (NaN unless sw.esr=99)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRAPH single solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution at layer iC of the game
iC=size(g,1);%end game

scrsz = get(0,'ScreenSize');
figure1 = figure('Name','Endgame Solution','Color',[1 1 1],'Position',[1 1 scrsz(3) scrsz(4)]);
graph.solution(g,iC,par,1, figure1);
if sw.alternate==1
    figure2 = figure('Name','Endgame Solution by turn','Color',[1 1 1],'Position',[1 1 scrsz(3) scrsz(4)]);
    graph.solution(g,iC,par,2, figure2);
end
% Edges at layer iC of the game
graph.edges(g,iC, par);

% Best responce functions
% ONLY works for simultanous move game - when br is returned
% parameters of interest:
%c1=(par.cmax-par.cmin)/2;
%c2=(par.cmax-par.cmin)/2;
%iC=size(g,1);%end game
%step=.001;
%graph.br(br,g,mp.eta,c1,c2,iC,step);

% Second order best responce functions
% (by roots)
% same arguments + firm index
%iF=1;
%graph.br2byroots(br,g,mp.eta, c1,c2,iC,iF,step)

% Simulation graphs
% NOT WORKING AT THE MOMENT




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOLVE the model for a ALL FEASIBLE ESR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sw.esr=99;

% If sw.esr is ESRstring (sw.esr=99) the loop is started over all possible 
% equilibrium selection rules which are generated in lexicographical order
% as sequences of digits (0 to MAXEQB-1) of length neqstr=nC*(nC+1)*(2*nC+1)/6
% where nC is number of grid points, and neqstr is essentially a sum of
% squares of natural numbers between 1 and nC
% Denote eqbstr one such sequence, let's call it ESR string
% ESR strings should be understood in the following way: the lowest=right-most digit 
% denotes the selected equilibrium at the top layer of the game (where iC=ic1=ic2=nC)
% Next four digits represent ESR in each point of the second layer of the game: 
% second lowest digit is interior point, next two digits are edges, and fourth lowest
% digit is the corner of the iC=nC-1 stage game.
% Next nine digits represent ESR on the third layer of the game (iC=nC-2) using the
% same block structure, namely lower 4 points correspond to interior, next 2 points
% represent one edge, 2 points another edge, and last point - the corner.
% Next 16 digits represent ESR on the fourth layer, and so on.
% 
% The loop over ESR strings works like this: first initial run is done using (0,..0)
% string (which is identical to selecting first available equilibrium in every point
% of the state space). This creates the game structure G which holds information on the
% number of stage equilibria in every point. The principle of indexing ESR strings
% ensures the following property: for each point in the state space (iC,ic1,ic2) 
% corresponding to digit K it holds that all points that may be dependent on the 
% equilibrium chosen at (iC,ic1,ic2) correspond to LOWER digits than K. Therefore 
% for any given ESR string after the initial run the game is recalculated only at 
% the state point corresponding to the highest non-zero digit, and up.

% Two more parameters must be set:
% Maximum number of DISTINCT outcomes from different ESR to be outputed (used for 
% allocation of memory). Disctinct values of ALL statistics.
sw.esrmax=100000;
% Index of the first lexicographical ESR string to be used
sw.esrstart=0; %index of the first eqstring

% run the solver
fprintf('Running the solver.. ');tic;
[bne, br, g, eqbstr]=leapfrog(par,mp,sw);
ts=toc;fprintf('done in %1.10fsec\n',ts);
% Output is:
%   per-period payoffs
%   coefficients from analytical solution (ONLY in simultanous move game)
%   game structure g
%   table with equilibrium statistics

% reshape the eqbstr output (delete NaN colums in eqbstrings and transpose)
eqbstr(:,isnan(eqbstr(1,:)))=[];
eqbstr=eqbstr';

% Equilibruim statistics in eqbstr are now the following:
% 1 column: lexicographical index of first ESR that gives outcome in this row
% 2 column: number of these outcomes (distinct pairs of value functions)
% 3 column: value of not investing of firm1
% 4 column: value of not investing of firm2 (when it's not its turn in alteranting move game)
% 5 column: whether equilibrium is pure in all state points
% 6 column: whether equilibrium is symmetric in all value functions
% 7 column: whether equilibrium includes a positive prob of investment of higher cost firm
%           in any state point (leapfrogging)
% These statistics can be edited in the outputiter() routine in the leapfrog.c 

isymmetric=eqbstr(eqbstr(:,6)==1,1); %save the index of symmetric ESR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRAPH collections of equilibria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: 
%   eqbstr - must contain statistics in columns
%   [s c]  - columns of statistics to be shown by size and color of the dots
%            when only size is given, B/W is drawn
%            when [] is given, standard size and B/W is drawn 
%   value of monopoly outcome, no line when [] is given

% calculate the monopoly outcome for current parameters
m_mp=mp;
m_mp.tpm=[1 0;1 0];
m_sw=sw;
m_sw.alternate=true;
m_sw.esr=5;
fprintf('Running the solver for monopoly game.. ');tic;
[a,a,mon,a]=leapfrog(par,mp,sw);
ts=toc;fprintf('done in %1.10fsec\n',ts);
% mon(par.nC).solution(1,7) is the value of not investing by firm 1

% build the graph in color
% 2--> use count for size and color
[data labels]=graph.EqbstrPlot(eqbstr,[2 2],mon(par.nC).solution(1,7));
title('Number and structure of equilibria');

% build the B/W graph
% 5--> use pure strategy indicator
[data labels]=graph.EqbstrPlot(eqbstr,5,mon(par.nC).solution(1,7));
title('Pure strategies');

% Output from graph.EqbstrPlot contains the statistics for the equilibria description
% table. labels contain some labels for the output, currently:
% 1 Total number of equilibria
% 2 Number of distinct pairs of value functions
% 3 Number of pure equilibria
% 4 Number of symmetric equilibria
% 5 Number of leapfrogging equilibria
out.total=data(1);
out.distinct=data(2);
out.pure=data(3);
out.symmetric=data(4);
out.leapfrogging=data(5);
out.maxValue_firm1=max(eqbstr(:,3));
out.maxValue_firm2=max(eqbstr(:,4));
out.monopoly=mon(par.nC).solution(1,7);
fprintf('Equilibria found:\n');
disp(out);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOLVE the model for a particular ESR string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sw.esr=99;

% To solve the model for a particular ESR string it is enough to set the 
% sw.esrstart and sw.esrmax parameters
% When sw.esrmax==1 the loop over ESR string is in special mode and only performs
% one iteration - to the next feasible ESR string after the one given.
% This is in contrast to sw.esrmax>1 when ALL equilibria will be computed and the
% number of them reported, even if sw.esrmax is too small to ouput all found 
% equilibria in the last output argument of leapfrog (eqbstr)

sw.esrmax=1;
sw.esrstart=5643; %index of the first eqstring
fprintf('Running the solver.. ');tic;
[a,a,g,eqbstr]=leapfrog(par,mp,sw);
ts=toc;fprintf('done in %1.10fsec\n',ts);
fprintf('Game structure g contains solution for the following ESR string:\n');
disp(lexistring(eqbstr(1),3)); % 3 must match MAXEQB value at compile time !!!

% Matlab function lexistring converts lexi-index to lexi-string and vise versa
% for a given MAXEQB (module)

% Saved earlier uique symmetric equilibrium is obtained from the string
disp(lexistring(isymmetric,3));

% To solve the model for a given ESR string
esr=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 0 0 0 2 0];%1+4+9+16 elements for nC=4
% run:
sw.esrstart=lexistring(esr,3);
fprintf('Running the solver.. ');tic;
[a,a,g,eqbstr]=leapfrog(par,mp,sw);
ts=toc;fprintf('done in %1.10fsec\n',ts);
fprintf('Game struncture g contains solution for the following ESR string:\n');
disp(lexistring(eqbstr(1),3)); % 3 must match MAXEQB value at compile time !!!

% It appears that there are issues with large lexi-indeces that are passed to C code
% as double floating point. Once even one digit is lost from the index it is not possible
% to reconstruct the string, and the solution may return NaN in eqbstr indicating there is 
% no feasible equilibrium past the index specified which could resulted from this 
% round off problem.


