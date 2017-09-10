%  Results.m: generates graphs for section 4 in Iskhakov, Rust and Schjerning (IER, 2017)
%         Fedor Iskhakov, Australian National University
%         John Rust, Georgetown University
%         Bertel Schjerning, University of Copenhagen

%% CLEAR MEMORY AND SWT SWITCHES
clc; clear; close all;

% to be modified by user
%paperdir='/Users/bertelschjerning/Dropbox/1PROJECTS/IRS/leapfrogging/_IOpaper/IER_figures/';
paperdir='results/';

%%Switches
compile=1;
fig5=0;
fig6=1
fig89=0;
loadresults=0;
saveresults=1;


%% COMPILE
tic
if compile
    fprintf('COMPILING... ');
    tic
    mex -largeArrayDims leapfrog.c -DMAXEQB=5 -DPRINTeqbloop=1 -DPRINTeqbstr=0
    tc=toc; fprintf('Compiled model in %1.10f (seconds)\n\n',tc);
end

% SET BECNHMARK PARAMETERS
setup;

%% FIGURE 5: Payoff maps and efficiency of simulateneous move game - with and without random technology
if fig5
    clear mpvec swvec parvec

    for i=1:2;
        mpvec(i)=mp;
        swvec(i)=sw;
        parvec(i)=par;

        swvec(i).alternate=false;
        parvec(i).nC=5; 
        mpvec(i).k1=8.3;
        mpvec(i).k2=1;

        swvec(i).esr=99; %equilibrium selection rule to be used: see setup.m or esr.c
        swvec(i).esrmax=30000000; % max number of district equilirbium payoffs
        swvec(i).esrstart=0;  %
        
    end
    
    i=1; % Deterministic technological progress
    tit1(i)={'Panel (a): Deterministic technological progress'}
    mpvec(i).onestep=1; 
    mpvec(i).c_tr=-1; 

    i=2; % Random technological progress
    tit1(i)={'Panel (b): Random technological progress'}
    mpvec(i).onestep=0; 
    mpvec(i).c_tr=-1; 

    
    for i=1:numel(mpvec);
        if swvec(i).alternate==false
            filenamefig5='fig5_sm_det.mat';
        else
            filenamefig5='fig5_sm_det.mat';
        end
        
        
        % recalculate depend parameters
        mpvec(i).nC=parvec(i).nC;
        parvec(i).pti=f_pti(mpvec(i), parvec(i));
        
        
        % SET MONOPOLY PARAMETERS
        mp_mon=mpvec(i);
        par_mon=parvec(i);
        sw_mon=swvec(i);
        mp_mon.tpm=[1 0; 1 0];
        sw_mon.alternate=true; % always  solve alternate move to obtain monopoly profits
        sw_mon.esr=1; %equil5ibrium selection rule to be used: see setup.m or esr.c
        sw_mon.esrmax=1; %run 1 eqstrings
        
        % SOLVE MONOPOLY MODEL
        [bne_mon, br_mon, g_mon, eqbstr_mon]=leapfrog(par_mon,mp_mon,sw_mon);       
        v10_mon=g_mon(par_mon.nC).solution(1,7);
        v20_mon=g_mon(par_mon.nC).solution(1,9);
        
        % SOLVE DUOPOLY MODEL AT PARAMETRS mpvec(i), parvec(i), swvec(i)
        if loadresults
            disp('Load results from duopoly game')
            load(filenamefig5)
            loadresults=1; % maintain loadresults=1 if results are loaded.
            saveresults=0; % do not overwrite loaded results;
        else
            disp('Solve dynamic duopoly game')
            [bne, br, g, eqbstr, moncmp]=leapfrog(parvec(i),mpvec(i),swvec(i),g_mon);
            eqbstr(:,isnan(eqbstr(1,:)))=[];
            eqbstr=eqbstr';
            moncmp=moncmp';
            %Efficiency
            eqbstr(:,8)=eqbstr(:,8)/(v10_mon+v20_mon);

            if saveresults
                save(filenamefig5);
            end
        end
                
        [a, b]=graph.EqbstrPlot(eqbstr,[2 8],(v10_mon+v20_mon),[],swvec(i),tit1{i});
        
        set(gcf, 'PaperPosition', [0 0 6 6]);
        saveas(gcf, sprintf('%sFig5_%d',paperdir,  i), 'fig')
                
        d=analyse(eqbstr);
    end
end



%% FIGURE 6: Payoff maps and efficiency of MPE in two specifications of the game
if fig6
    clear mpvec swvec parvec
    for i=1:2;
        mpvec(i)=mp;
        swvec(i)=sw;
        parvec(i)=par;

        parvec(i).nC=5;
        mpvec(i).c_tr=1; 
        mpvec(i).onestep=1; 

        mpvec(i).k1=5;
        mpvec(i).k2=0;
        mpvec(i).tpm=[0 1; 1 0];


        swvec(i).esr=99; %equilibrium selection rule to be used: see setup.m or esr.c
        swvec(i).esrmax=30000000; % max number of district equilibrium payoffs 
        swvec(i).esrstart=0;  %
        
    end
    
    i=1; %% alternating move
    tit1(i)={'Panel (a): Alternating move'}
    tit2(i)={'Panel (b): Alternating move'}
    swvec(i).alternate=true;
        
    i=2; 
    tit1(i)={'Panel (c): Simultaneous move'}
    tit2(i)={'Panel (d): Simultaneous move'}
    swvec(i).alternate=false;
    
    for i=1:numel(mpvec);
        if swvec(i).alternate==false
            filenamefig6='fig6_sm.mat';
        else
            filenamefig6='fig6_am.mat';
        end
        
        
        % recalculate depend parameters
        mpvec(i).nC=parvec(i).nC;
        parvec(i).pti=f_pti(mpvec(i), parvec(i));
        
        
        % SET MONOPOLY PARAMETERS
        mp_mon=mpvec(i);
        par_mon=parvec(i);
        sw_mon=swvec(i);
        mp_mon.tpm=[1 0; 1 0];
        sw_mon.alternate=true; % always  solve alternate move to obtain monopoly profits
        sw_mon.esr=1; %equilibrium selection rule to be used: see setup.m or esr.c
        sw_mon.esrmax=1; %run N feasible eqstrings (set equal to 192736405 for n=5)
        
        % SOLVE MONOPOLY MODEL
        [bne_mon, br_mon, g_mon, eqbstr_mon]=leapfrog(par_mon,mp_mon,sw_mon);       
        v10_mon=g_mon(par_mon.nC).solution(1,7);
        v20_mon=g_mon(par_mon.nC).solution(1,9);
        
        % SOLVE DUOPOLY MODEL AT PARAMETRS mpvec(i), parvec(i), swvec(i)
        if loadresults
            disp('Load results from duopoly game')
            load(filenamefig6)
            loadresults=1; % maintain loadresults=1 if results are loaded.
            saveresults=0; % do not overwrite loaded results;
        else
            disp('Solve dynamic duopoly game')
            [bne, br, g, eqbstr, moncmp]=leapfrog(parvec(i),mpvec(i),swvec(i),g_mon);
            eqbstr(:,isnan(eqbstr(1,:)))=[];
            eqbstr=eqbstr';
            moncmp=moncmp';
            %Efficiency
            eqbstr(:,8)=eqbstr(:,8)/(v10_mon+v20_mon);

            if saveresults
                save(filenamefig6);
            end
        end
        
        
        [a, b]=graph.EqbstrPlot(eqbstr,[2 8],(v10_mon+v20_mon),[],swvec(i),tit1{i});
        
        set(gcf, 'PaperPosition', [0 0 6 6]);
        saveas(gcf, sprintf('%sTri_%d',paperdir,  i), 'fig')
                
        d=analyse(eqbstr);
        
        dummystr = '[((eqbstr(:,3)~=0).*(eqbstr(:,4)~=0)) ((eqbstr(:,3)==0) | (eqbstr(:,4)==0))]';
        graphstring='graph.EqbstrCDFPlot(eqbstr,[0], dummystr, lbl ,mp,sw,tit2{i},0.5, mfig)';
        %  graph.EqbstrCDFPlot(eqbstr,[0 -5 6 7], dummystr, {'vi ne 0','vi=0'} ,[],sw,tit2{i},1);
        graph.EqbstrCDFPlot(eqbstr,[0 7 -5 ], [], [] ,[],swvec(i),tit2{i},1);
        %   title(tit2{i});
        set(gcf, 'PaperPosition', [0 0 6 6]);
        saveas(gcf, sprintf('%sCDF_%d',paperdir,  i), 'fig')
    end
end

%% FIGURE 8 and 9
if fig89
    clear mpvec, swvec;
    mp.k1=2;
    for i=1:7;
        mpvec(i)=mp;
        swvec(i)=sw;
    end
    
    % Adjust mp structure - see setup for description
    % FIGURE 8
    i=1; mpvec(i).k1=2;  mpvec(i).nC=100; swvec(i).esr=99;		    tit(i)={'Panel (a): Preemption and rent-dissipation'}
    i=2; mpvec(i).k1=2;  mpvec(i).nC=25;  swvec(i).esr=99;		    tit(i)={'Panel (b): Underinvestment'}
    i=3; mpvec(i).k1=.5; mpvec(i).nC=25;  swvec(i).esr=99;		    tit(i)={'Panel (c): Leap-frogging'}
    
    
    % FIGURE 9
    i=4; mpvec(i).tpm=[0.2 0.8; 0.8 .2];  							tit(i)={'Panel (a): Random alternating moves'}
    i=5; mpvec(i).c_tr=1;  											tit(i)={'Panel (b): Non-monotonic tech. progress'}
    i=6; mpvec(i).c_tr=1; mpvec(i).onestep=0; 						tit(i)={'Panel (c): Non-monotonic multistep tech. progress'}
    i=7; mpvec(i)=mp; swvec(i).esr=7; swvec(i).alternate=false;		tit(i)={'Panel (d): Simultaneous move'}
    
    %!rm randstream.mat
    
    %% FIGURE 6 and 7
    for i=1:numel(mpvec);
        % for i=1:2;
        % recalculate depend parameters
        
        par.nC=mpvec(i).nC;  % add nC to mp
        capT=floor(1.5*par.nC);
        mpvec(i).dt=25/par.nC;            % 5 is top cost
        mpvec(i).df=exp(-0.05*mp.dt);
        par.pti=f_pti(mpvec(i), par);
        
        %% SET MONOPOLY PARAMETERS
        mp_mon=mpvec(i);
        par_mon=par;
        sw_mon=swvec(i);
        mp_mon.tpm=[1 0; 1 0];
        sw_mon.alternate=true; % always  solve alternate move to obtain monopoly profits
        sw_mon.esr=1; %equil5ibrium selection rule to be used: see setup.m or esr.c
        sw_mon.esrmax=1; % max number of distringt equilirbia (set equal to 16510 for sm with n=5)
        
        %% SOLVE MONOPOLY MODEL
        disp('solve monopoly model')
        [bne_mon, br_mon, g_mon, eqbstr_mon]=leapfrog(par_mon,mp_mon,sw_mon);
        
        v10_mon=g_mon(par.nC).solution(1,7);
        v20_mon=g_mon(par.nC).solution(1,9);
        
        %% SOLVE MODEL AT PARAMETRS mpvec(i), par, swvec(i)
        [bne, br, g, eqbstr]=leapfrog(par,mpvec(i),swvec(i));
        
        %% Simulare sequences
        sp=simul.setup(par);         % Run setup for simulation module
        sp.T=capT;
        s=simul.sequence(g,sp, mpvec(i),par,swvec(i).alternate);   % Sumulate sequences
        s_mon=simul.sequence(g_mon,sp, mp_mon,par,sw_mon.alternate);   % Sumulate sequences
        
        
        
        if i<=7;
            s.cmon=s_mon.c1;
            graph.CostSequence(s, mpvec(i), tit{i}, '', s_mon); % plot sequences realized costs
        else
            graph.CostSequence(s, mpvec(i), tit{i}, ''); % plot sequences realized costs
        end
        set(gcf, 'PaperPosition', [0 0 6 5]);
        saveas(gcf, sprintf('%sSeq_%d',paperdir,  i), 'fig')
        saveas(gcf, sprintf('%sSeq_%d',paperdir,  i), 'pdf')
        saveas(gcf, sprintf('%sSeq_%d.eps',paperdir,  i), 'psc2')
    end
end



