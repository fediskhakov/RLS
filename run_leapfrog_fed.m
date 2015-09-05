%         Fedor Iskhakov, University Technology Sidney
%         John Rust, University of Maryland
%         Bertel Schjerning, University of Copenhagen
%         Sidney, March 2012

close all;

%% COMPILE 
%if false
if true
    fprintf('Compiling... '); 
    tic
    mex -largeArrayDims leapfrog.c -DMAXEQB=3 -DPRINTeqbloop=1 -DPRINTeqbstr=0
    tc=toc; fprintf('Compiled model in %1.10f (seconds)\n\n',tc);
end

%% PARAMETERS
% load defaults
setup;
% > time differential
%mp.dt=1;
% > investment cost
mp.k1=5.2;
mp.k2=0;
% > eta
mp.eta=0;
%mp.sigma=10;
% > technological progress
mp.c_tr=-.75;
mp.onestep=1;  
mp.beta_a=0.8; 
mp.beta_b=0.4;           
% > transition probability for turns (standard indexing: from in rows, to in colums)
%mp.tpm=[.1 .9;
%        .9 .1]; 
mp.tpm=[0 1;
        1 0]; 
% > number of grid points
par.nC=150;

%LOAD parameters from file
%load('../runs/SM-nc-eta-pr-step/results_run00061.mat');
%desc='(1) DTP';
%mp=mp_run00061;
%par=par_run00061;
%sw=sw_run00061;
%clear *_run00061;
%load('../runs/SM-nc-eta-pr-step/results_run00062.mat');
%desc='(2) MSTP';
%mp=mp_run00062;
%par=par_run00062;
%sw=sw_run00062;
%clear *_run00062;
desc='';

% recalculate depend parameters
[par mp]=f_update_params(par,mp);

%MONOPOLY SOLUTION
mpm=mp;
mpm.tpm=[1 0; 1 0];
swm.alternate=true;
swm.analytical=true;
swm.esr=5;
% swm.esr=99;
% swm.esrstart=0;
% swm.esrmax=1;
swm.esrmax=1;
swm.esrstart=0;
[a,b,gmon,c]=leapfrog(par,mpm,swm);
mon=gmon(par.nC).solution(1,7)+gmon(par.nC).solution(1,9);
sp=simul.setup(par); %run setup for simulation module
sp.T=floor(1.5*par.maxit);
simmon=simul.sequence(gmon,sp,mpm,par,1);   % Sumulate sequences
clear a b c swm sp;

%% MODEL SWITCHES
sw.alternate=true; %alternate move or simultanious move game
sw.analytical=true; %analytical solution of eta>0 simultanous move game
sw.esr=99; %equilibrium selection rule to be used: see setup.m or esr.c
%sw.esr=5; %ESR1
%sw.esr=3; %HighCostInvestMix
sw.esrmax=100000; %run N feasible eqstrings
%sw.esrmax=1; %run N feasible eqstrings
sw.esrstart=0; %index of the first eqstring

 % sw.esrstart=243;
 % sw.esrmax=1;

%sw.esrstart=[0 0 0 0 0 0 0 1 2 2 0 1 2 0 0 1 0 0 0 0 0 1 2 0 1 0 0 0 1 0]; %leapfrogging
%sw.esrstart=[0 0 0 0 0 0 0 1 0 0 2 1 0 2 2 1 0 0 0 0 0 1 0 2 1 0 0 0 1 0]; %anti-leapfrogging
%sw.esrstart=[0 0 0 0 0 0 0 1 0 0 2 1 0 2 2 1 0 0 0 0 0 1 0 0 1 0 0 0 0 0]; %anti-leapfrogging as ran by esr=6
%sw.esrmax=1;


%% RUN LRS LOOP ONCE
%if false
if true
    fprintf('Parameters:\n')
    disp(par);
    disp(mp);
    disp(sw);
    fprintf('Transition probability for m:\n')
    disp(mp.tpm);
    
    tic; 
    [bne, br, g, eqbstr, moncmp]=leapfrog(par,mp,sw,gmon); %monopoly g structure for underinvestment checks
    ts=toc;
    fprintf('Solved model in %1.10f (seconds)\n',ts);
    eqbstr(:,isnan(eqbstr(1,:)))=[];
    eqbstr=eqbstr';
    moncmp=moncmp';

    %Efficiency
    %eqbstr(:,8)=1-eqbstr(:,8)/mon;
    eqbstr(:,8)=eqbstr(:,8)/mon;
    if sw.alternate
        d=analyse(eqbstr,'alt');
        d=analyse(eqbstr,'notab','alt');
    else
        d=analyse(eqbstr);
        d=analyse(eqbstr,'notab');
    end
    eqbsum=graph.EqbstrPlot (eqbstr,[2 2],mon,mp,sw,'Pay-off map for all found equilibria');

    %for the figures in the paper:
    %skip long title by passing [] in place of mp
    graph.EqbstrPlot(eqbstr,[8 8],mon,mp,sw,desc);
    %graph.EqbstrPlot(eqbstr,2,mon,[],sw,desc);
    %graph.EqbstrPlot(eqbstr,[5 2],mon,[],sw,['(2) ' desc ': pure equilibria']);
    %graph.EqbstrPlot(eqbstr,[6 2],mon,[],sw,['(3) ' desc ': symmetric equilibria']);


    if 1
        %SIMULATE
        sp=simul.setup(par);         % Run setup for simulation module
        sp.T=floor(1.5*par.maxit);
        s=simul.sequence(g,sp,mp,par,sw.alternate);   % Sumulate sequences
        graph.CostSequence(s, mp, 'Duopoly equilibrium realization'); % plot sequences realized costs
        % graph.ProfitSequence(s, mp, 'Duopoly profits'); % plot sequences realized costs
        graph.CostSequence(simmon, mpm, 'Monopoly equilibrium realization'); % plot sequences realized costs
    end
end

%% RUN MODLE ONCE and save results to disk
if false
%if true
    [bne, br, g, eqbstr]=RunLeap(mp,par,sw,'-DMAXEQB=3 -DPRINTeqbloop=0 -DPRINTeqbstr=0');
end

%% RUN MovieMaker
if false
%if true
    dr ='../../../movies/';
    %dr='';
    %moviemaker (par,mp,sw,'mp.k1',[0 100],[dr 'eqb.k1.n5']);
    %moviemaker (par,mp,sw,'mp.tpm(1,1)',[1 0],[dr 'eqb.tpm.n5']);
    moviemaker (par,mp,sw,'mp.eta',[0.5 0],[dr 'eqb.eta.n5.half2zero'],[1 12]);
    %moviemaker (par,mp,sw,'mp.df',[.8 .99999],[dr 'eqb.df.n5']);

    %screen play
    %moviemaker (par,mp,sw,'mp.eta',[.5 0],'',[1 12]);
    %moviemaker (par,mp,sw,'mp.beta_a',[0.001 .5],'');
    %moviemaker (par,mp,sw,'mp.c_tr',[0.001 .1],'');

end

clear tc ts;


