%         Fedor Iskhakov, University Technology Sidney
%         John Rust, University of Maryland
%         Bertel Schjerning, University of Copenhagen
%         Sidney, March 2012

% This script creates 3 graphs of best response functions for the Appendix
% SET how many stage equilibria to search for in mex statement PRINTnrstagempe=5 (line 18)

%% CLEAR MEMORY AND SET SWITCHES
clc
clear;
close all;

%% COMPILE 
if 1
    tic
    fprintf('Compiling... '); 
    tic
    mex -largeArrayDims leapfrog.c -DMAXEQB=5 -DPRINTeqbloop=0 -DPRINTeqbstr=0 -DPRINTnrstagempe=5
    tc=toc; 
    fprintf('Compiled model in %1.10f (seconds)\n\n',tc);
end

%% SET BECNHMARK PARAMETERS 
setup;
% Adjust mp structure - see setup for description
% mp.k1=10; mp.k2=1.0; mp.c_tr=1;  % ASSYMMETRIC
% mp.k1=1; mp.k2=10; mp.c_tr=0.05; % SYMMETRIC
mp.eta=0;
mp.tpm=[0 1;
        1 0]; 
mp.c_tr=-.5;
mp.c_tr=0.05
mp.beta_a=1.8; 
mp.beta_b=0.4;    

% Adjustments to par structure - see setup for description
par.nC=4;
par.pti=f_pti(mp, par);

% recalculate depend parameters
[par mp]=f_update_params(par,mp);

% Adjustments to sw structure - see setup for description
% IDEA: run whole RLS to find ESR where we have desirable number of stage equilibria (see DPRINTnrstagempe switch in compile)
sw.alternate=false; %alternate move or simultanious move game
sw.esr=99; %equilibrium selection rule to be used: see setup.m or esr.c
sw.esrmax=1000000; %run N feasible eqstrings (set equal to 192736405 for n=5)
sw.esrstart=0; %index of the first eqstring

fprintf('Parameters are:\n');
mp
tpm=mp.tpm
par


if 0
% Run this block to output ESR string that have 5 equilibria stages for the next part if needed

    %% PART 1: full RLS with output of ESR strings with particular number of stage equilibria given in -DPRINTnrstagempe=5

    fprintf('Solving model... \n')
    tic;
    [bne, br, g, eqbstr]=leapfrog(par,mp,sw);
    ts=toc;
    fprintf('Solved model in %1.10f (seconds)\n',ts);
    eqbstr(:,isnan(eqbstr(1,:)))=[];
    eqbstr=eqbstr';

    %% Plot of all equilibria (without monopoly values)
    scrsz = get(0,'ScreenSize');
    fig1 = figure('Name','All equilibria','Color',[1 1 1],'Position',[1 1 scrsz(3) scrsz(4)]);
    [a, b]=graph.EqbstrPlot(eqbstr,[2 5],[],mp, sw,'All eqb',fig1);
    % 1: lex index of
    % 2: number of repetions (sum of col2 numer of equilibria in total)
    % 3: v10
    % 4: v20 (x20 if alternating move)
    % 5: stat1 pure strategy dynamic equilibrium (0/1)
    % 6: stat2 symetric dynamic aquilibrium (0/1)
    % 7: stat3 leapfrogging eqb (high cost follower has poitsitive investment probability) (0/1)



else %run full RLS or only one ESR string

    %% PART 2: running for particular ESR and making plots

    % Solve the model with the ESR that has 5 equilibria stages
    sw.esrstart=[0 0 0 0 0 0 0 2 2 0 2 2 1 0 2 1 0 0 0 0 0 2 0 2 2 0 0 0 4 0];
    %Other parameters
    sw.esrmax=1;%only run solver onces
    %run solver
    [bne,br,g,eqbstr]=leapfrog(par,mp,sw);

    %draw plot for the following (ic ic1 ic2):
    doplots=           {[1 1 1]};     % 1 stage game MPE 
    doplots={doplots{:},[1 3 2]};     % 3 (2 pure + 1 mixed)
    doplots={doplots{:},[2 4 4]};     % 3 (1 pure + 2 mixed)
    doplots={doplots{:},[3 4 4]};     % 5 
    %make codes for quick checks
    for i=1:numel(doplots)
        doplotcodes(i)=doplots{i}*[100 10 1]';
    end

    %make doplots={} to do all plots
    doplotcodes={};

    % Report number of stage equilibria on all stages
    fprintf('\nNumber of stage equilibria:\n');
    fprintf('%s\n',repmat('-',1,44));
    fprintf('%3s %3s %3s %10s %6s %6s %6s\n','ic','ic1','ic2','type','pure','mixed','total');
    fprintf('%s\n',repmat('-',1,44));
    for ic=1:par.nC
    % for ic=2        
        gq=g(ic).solution(:,[1 2 11 12]); %take ic1 ic2 p1 p2 from the solution structure
        for ic1=ic:par.nC
            for ic2=ic:par.nC
                if ic1==ic & ic2==ic
                    label='corner';
                elseif ic1==ic | ic2==ic
                    label='edge';
                else
                    label='interior';
                end
                indx=find(gq(:,1)==ic1-1 & gq(:,2)==ic2-1);
                numeqb=numel(indx);
                numpure=0;
                nummix=0;
                for i=reshape(indx,1,[])
                    if (abs(gq(i,3))<eps || abs(gq(i,3)-1)<eps) && (abs(gq(i,4))<eps || abs(gq(i,4)-1)<eps)
                        numpure=numpure+1;
                    else
                        nummix=nummix+1;
                    end
                end
                if numeqb==5 || (~isempty(doplotcodes) && ismember([ic ic1 ic2]*[100 10 1]',doplotcodes))
                    mark=' <<<';
                else
                    mark='';
                end
                %entry into the table
                fprintf('%3d %3d %3d %10s %6d %6d %6d%s\n',ic,ic1,ic2,label,numpure,nummix,numeqb,mark);
                %plot best responses (for all state points)
                if isempty(doplotcodes) || ismember([ic ic1 ic2]*[100 10 1]',doplotcodes)
                    if mp.eta<eps
                        graph.br0eta(br,g,ic1,ic2,ic);
                    else
                        c_grid=linspace(par.cmin,par.cmax,par.nC);
                        spacing=0.001;
                        c1=c_grid(ic1);
                        c2=c_grid(ic2);
                        graph.br(br,g,mp.eta,c1,c2,ic,spacing);
                    end
                end
            end
        end
        fprintf('%s\n',repmat('-',1,44));
    end
    fprintf('\n');
end
