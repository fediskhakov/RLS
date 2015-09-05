%         Fedor Iskhakov, University Technology Sydney
%         John Rust, Georgetowb University
%         Bertel Schjerning, University of Copenhagen
%         Sliver Sping, August 2012

%% CLEAR MEMORY AND SWT SWITCHES
clc; clear; close all;

% Switches (Local to this script)
compile=1;
printparam=0; %% switch of to make output more compact
capT=440;   %% max number of time periods in time simulations

%% COMPILE
tic
if compile
    fprintf('COMPILING... ');
    tic
    %     mex -largeArrayDims leapfrog.c -DMAXEQB=3 -DPRINTeqbloop=1 -DPRINTeqbstr=0
    mex -largeArrayDims leapfrog.c -DMAXEQB=3
    tc=toc; fprintf('Compiled model in %1.10f (seconds)\n\n',tc);
end

%% SET BECNHMARK PARAMETERS
setup;

% Adjust mp structure - see setup for description
mp.eta=0;
mp.tpm=[0 1;
        1 0];
%mp.tpm=[0.5 0.5;
%        0.5 0.5];

mp.k1=10;
mp.k2=0;
par.nC=150;                 % 4
mp.dt=80/par.nC;           % 5
mp.df=exp(-0.05*mp.dt);
mp.c_tr=-.5;
mp.onestep=1;

% recalculate depend parameters
par.pti=f_pti(mp, par);
mp.nC=par.nC;  % add nC to mp

% Adjustments to sw structure - see setup for description
sw.alternate=true; %alternate move or simultanious move game
sw.esr=5; %equilibrium selection rule to be used: see setup.m or esr.c
sw.esrmax=1000000; %run N feasible eqstrings (set equal to 192736405 for n=5)
sw.esrstart=0; %index of the first eqstring

% Run particular eqbstring (Uncommen next thee lines to run single ESR)
% sw.esr=99; %equilibrium selection rule to be used: see setup.m or esr.c
% sw.esrmax=1; %run N feasible eqstrings (set equal to 192736405 for n=5)
% sw.esrstart=2324522934;  % (number stored second column of eqbstr)

%% SET MONOPOLY PARAMETERS
mp_mon=mp;
par_mon=par;
sw_mon=sw;

% Adjustments to par structure - see setup for description
mp_mon.tpm=[1 0; 1 0];

% Adjustments to sw structure - see setup for description
sw_mon.alternate=true; % always  solve alternate move to obtain monopoly profits
sw_mon.esr=1; %equil5ibrium selection rule to be used: see setup.m or esr.c
sw_mon.esrmax=1; %run N feasible eqstrings (set equal to 192736405 for n=5)

%% SOLVE MONOPOLY MODEL
fprintf('MONOPOLY MODEL\n');
% set parameters to benchmark parameters
if printparam
    fprintf('Parameters are:\n');
    mp_mon
    tpm=mp_mon.tpm
    par_mon
    sw_mon
end

%tic; [bne_mon, br_mon, g_mon, eqbstr_mon]=leapfrog(par_mon,mp_mon,sw_mon); ts=toc;
%fprintf('Solve model in %1.10f (seconds)\n\n\n\n',ts);

%mon=g_mon(par.nC).solution(1,7)+g_mon(par.nC).solution(1,9);

%% SOLVE MODEL AT PARAMETRS mp, par, sw
fprintf('BENCHMARK MODEL\n');
if printparam
    fprintf('Parameters are:\n');
    mp
    tpm=mp.tpm
    par
    sw.analytical=0
end

tic;
[bne, br, g, eqbstr]=leapfrog(par,mp,sw);
ts=toc;
fprintf('Solved model in %1.10f (seconds)\n',ts);
fprintf('doing maitnance to solution results\n\n');
eqbstr(:,isnan(eqbstr(1,:)))=[];
eqbstr=eqbstr';


% eqbstr:
% 1: lex index of
% 2: number of repetions (sum of col2 numer of equilibria in total)
% 3: v10
% 4: v20 (x20 if alternating move)
% 5: stat1 pure strategy dynamic equilibrium (0/1)
% 6: stat2 symetric dynamic aquilibrium (0/1)
% 7: stat3 leapfrogging eqb (high cost follower has poitsitive investment probability) (0/1)

%% PRINT MONOPOLY AND DUOPOLY PROFITS
%v10_mon=g_mon(par.nC).solution(1,7);
%v20_mon=g_mon(par.nC).solution(1,9);

v10_duo=g(par.nC).solution(1,7);
v20_duo=g(par.nC).solution(1,9);
if sw.alternate;
    v20_duo=g(par.nC).solution(1,18);
end;
fprintf('-------------------------------------------------------------------------------------\n');
%fprintf('Monopoly profits: %f\n', v10_mon+v20_mon);
fprintf('\nDuopoly  profits\n');
fprintf(' Firm 1: %f\n', v10_duo);
fprintf(' Firm 2: %f\n', v20_duo);
fprintf(' Total : %f\n', v10_duo+v20_duo);
fprintf('-------------------------------------------------------------------------------------\n');
fprintf('\n\n');



%% SIMULATE COST SEQUENCES
% Simulate and plot realized sequences
sp=simul.setup(par);         % Run setup for simulation module
sp.T=capT;
s=simul.sequence(g,sp, mp,par,0);   % Sumulate sequences
graph.CostSequence(s, mp, 'Duopoly equilibrium realization'); % plot sequences realized costs
graph.ProfitSequence(s, mp, 'Duopoly profits'); % plot sequences realized costs

%s_mon=simul.sequence(g_mon,sp, mp_mon,par_mon,1);   % Sumulate sequences
%graph.CostSequence(s_mon, mp, 'Monopoly equilibrium realization'); % plot sequences realized costs
%graph.ProfitSequence(s_mon, mp_mon, 'Monpoly profits'); % plot sequences realized costs

%% All equilirbia
if (mp.nC<=5 || (sw.alternate & mp.nC<=20));
    fprintf('SOLVE BENCHMARK MODEL - for all equilribria\n');
    sw.esr=99; %equilibrium selection rule to be used: see setup.m or esr.c
    sw.esrmax=100000; %run N feasible eqstrings (set equal to 192736405 for n=5)
    sw.esrstart=0;  %
    
    tic;
    [bne, br, g, eqbstr]=leapfrog(par,mp,sw);
    ts=toc;
    fprintf('Solved model in %1.10f (seconds)\n',ts);
    fprintf('doing maitnance to solution results\n');
    eqbstr(:,1:min(10,size(eqbstr,2)));
    eqbstr(:,isnan(eqbstr(1,:)))=[];
    eqbstr=eqbstr';
    
    [a, b]=graph.EqbstrPlot(eqbstr,[2 5],[v10_mon+v20_mon],mp, sw,'All eqb');

    fprintf('\n\n');
    fprintf('Equilibria summary, all equilirbia\n');
    fprintf('-------------------------------------------------------------------------------------\n');
    for j=1:5; 
        fprintf('%-30s %15d\n', cell2mat(b(j)), a(j)); 
    end
    fprintf('-------------------------------------------------------------------------------------\n');
    fprintf('\n\n');
    
    %equilibrium strings where firm 1 and 2 achieve monopoly profits
    eqbstr_mon_f1=eqbstr(find(eqbstr(:,3)==v10_mon),:); % firm 1
    eqbstr_mon_f2=eqbstr(find(eqbstr(:,4)==v10_mon),:); % firm 2
    
    fprintf('Monopoly equilirbira\n');
    fprintf('      Lex index   # of ESR        v10    v20/x20       pure     syemmetric  leapfrog\n');
    fprintf('-------------------------------------------------------------------------------------\n');
    fprintf('%15d %10d %10g %10g %10d %10d %10d \n', eqbstr_mon_f1(:,1:7)')
    fprintf('%15d %10d %10g %10g %10d %10d %10d \n', eqbstr_mon_f2(:,1:7)')
    fprintf('\n');
    fprintf('Number of monopoly equilirbia   : %f\n', sum(eqbstr_mon_f1(:,2)+eqbstr_mon_f2(:,2)));
    fprintf('Number of equilirbia            : %f\n', sum(eqbstr(:,2)));
    fprintf('Monopoly equilibria (share of all equilirbia) : %f\n', sum(eqbstr_mon_f1(:,2)+eqbstr_mon_f2(:,2))/sum(eqbstr(:,2)));
    fprintf('-------------------------------------------------------------------------------------\n');
    
end;
