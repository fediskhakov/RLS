%  Results.m: generates graphs for section 4 in Iskhakov, Rust and Schjerning (2012)
%         Fedor Iskhakov, University of New South Wales
%         John Rust, University of Maryland
%         Bertel Schjerning, University of Copenhagen

%% CLEAR MEMORY AND SWT SWITCHES
clc; clear; close all;


%%Switches
compile=0;

%% COMPILE 
tic
if compile
fprintf('COMPILING... '); 
tic
%     mex -largeArrayDims$ leapfrog.c -DMAXEQB=3 -DPRINTeqbloop=1 -DPRINTeqbstr=0
mex -largeArrayDims leapfrog.c -DMAXEQB=3 -DPRINTeqbloop=1 -DPRINTeqbstr=0
tc=toc; fprintf('Compiled model in %1.10f (seconds)\n\n',tc);
end

% SET BECNHMARK PARAMETERS
setup;

% Adjustments to sw structure - see setup for description
sw.alternate=true; %alternate move or simultanious move game
sw.esr=5; %equilibrium selection rule to be used: see setup.m or esr.c %ESR:always play first equilibrium
mp.c_tr=-1;   % pc=dt*c_tr*c/(1.0+dt.c_tr*c);
mp.c_tr=1000;
mp.onestep=1;
mp.tpm=[0 1; 1 0];
mp.k1=5; 
mp.k2=0; 
par.nC=5;                 
mp.nC=par.nC;  % add nC to mp


%%%%%%%%%%%%%%%%%%%
%% payoff map
%%%%%%%%%%%%%%%%%%%


% recalculate depend parameters
mp.dt=5/par.nC;            % 5 is top cost
mp.df=exp(-0.05*mp.dt);
par.pti=f_pti(mp, par);

sw.esr=99; %equilibrium selection rule to be used: see setup.m or esr.c
sw.esrmax=100000; %run N feasible eqstrings (set equal to 192736405 for n=5)
sw.esrstart=0;  %

%% SET MONOPOLY PARAMETERS
mp_mon=mp;
par_mon=par;
sw_mon=sw;
mp_mon.tpm=[1 0; 1 0];
sw_mon.alternate=true; % always  solve alternate move to obtain monopoly profits
sw_mon.esr=1; %equil5ibrium selection rule to be used: see setup.m or esr.c
sw_mon.esrmax=1; %run N feasible eqstrings (set equal to 192736405 for n=5)

%% SOLVE MONOPOLY MODEL
[bne_mon, br_mon, g_mon, eqbstr_mon]=leapfrog(par_mon,mp_mon,sw_mon); 

v10_mon=g_mon(par.nC).solution(1,7);
v20_mon=g_mon(par.nC).solution(1,9);

%% SOLVE MODEL AT PARAMETRS mpvec(i), par, swvec(i)
[bne, br, g, eqbstr, moncmp]=leapfrog(par,mp,sw,g_mon);

eqbstr(:,isnan(eqbstr(1,:)))=[];
eqbstr=eqbstr';
moncmp=moncmp';

%Efficiency
%eqbstr(:,8)=eqbstr(:,8)/(v10_mon+v20_mon);
eqbstr(:,8)=eqbstr(:,8)/(v10_mon+v20_mon);

[a, b]=graph.EqbstrPlot(eqbstr,[2 8],(v10_mon+v20_mon),[],sw,'Payoff map alternating move game'); 

%%%%%%%%%%%%%%%%%%%
%% Movie 1
%%%%%%%%%%%%%%%%%%%
if false
dr=sprintf('%s/movies', pwd)
moviemaker (par,mp,sw,'mp.c_tr',[50 0],[dr 'eqb.c_tr_last50'],[]);
end

%%%%%%%%%%%%%%%%%%%
%% Movie 1
%%%%%%%%%%%%%%%%%%%
if true
dr=sprintf('%s/movies', pwd)
mp.c_tr=1;
moviemaker (par,mp,sw,'mp.eta',[1.5 0],[dr 'eqb.eta_1.5_to_0'],[]);
drawnow;
end

close all;

%%%%%%%%%%%%%%%%%%%
%% EQB correspondance v10 vs pi(c|c)
%%%%%%%%%%%%%%%%%%%
figure(1)
set(0,'DefaultAxesColorOrder',[0 0 0])	
for picc=[0.0000000:0.001:1];

    c_tr=(1/picc-1)/5;

	% recalculate depend parameters
	mp.c_tr=c_tr;
	par.pti=f_pti(mp, par);

	mp_mon.c_tr=c_tr;
	par_mon.pti=f_pti(mp_mon, par);

	%% SOLVE MONOPOLY MODEL
	[bne_mon, br_mon, g_mon, eqbstr_mon]=leapfrog(par_mon,mp_mon,sw_mon); 

	v10_mon=g_mon(par.nC).solution(1,7);
	v20_mon=g_mon(par.nC).solution(1,9);

	%% SOLVE MODEL AT PARAMETRS mpvec(i), par, swvec(i)
	[bne, br, g, eqbstr, moncmp]=leapfrog(par,mp,sw,g_mon);

	eqbstr(:,isnan(eqbstr(1,:)))=[];
	eqbstr=eqbstr';
	moncmp=moncmp';

	%Efficiency
	eqbstr(:,8)=eqbstr(:,8)/(v10_mon+v20_mon);

	v10=sort(eqbstr(:,3));
	v20=eqbstr(:,4);
    plot(picc,v10);
	hold on;
end;

    
    xlabel('\pi(c|c)')
	ylabel('vN1')
	title('Equilibrium correspondance (alternating move game), vN1 at (5,5,5) apex')
	set(gcf, 'PaperPosition', [0 0 6 5]);
    saveas(gcf, 'eqbcorrespondance_v10_pcc', 'fig')
    saveas(gcf, 'eqbcorrespondance_v10_pcc.eps', 'psc2')
disp('GOT HERE')
%%%%%%%%%%%%%%%%%%%
%% EQB correspondance v10 vs c_tr
%%%%%%%%%%%%%%%%%%%
figure(2)
set(0,'DefaultAxesColorOrder',[0 0 0])	
for eta=[1.5:-0.001:0];

    c_tr=(1/picc-1)/5;
    mp.c_tr=1; 

	% recalculate depend parameters
	mp.eta=eta;
	par.pti=f_pti(mp, par);

	mp_mon.c_tr=c_tr;
	par_mon.pti=f_pti(mp_mon, par);

	%% SOLVE MONOPOLY MODEL
	[bne_mon, br_mon, g_mon, eqbstr_mon]=leapfrog(par_mon,mp_mon,sw_mon); 

	v10_mon=g_mon(par.nC).solution(1,7);
	v20_mon=g_mon(par.nC).solution(1,9);

	%% SOLVE MODEL AT PARAMETRS mpvec(i), par, swvec(i)
	[bne, br, g, eqbstr, moncmp]=leapfrog(par,mp,sw,g_mon);

	eqbstr(:,isnan(eqbstr(1,:)))=[];
	eqbstr=eqbstr';
	moncmp=moncmp';

	%Efficiency
	eqbstr(:,8)=eqbstr(:,8)/(v10_mon+v20_mon);

	v10=sort(eqbstr(:,3));
	v20=eqbstr(:,4);
    plot(eta,v10);
	hold on;
end;

	xlabel('Investment shock scaling parameter, \eta')
	ylabel('vN1')
	title('Equilibrium correspondance (alternating move game), vN1 at (5,5,5) apex')
	set(gcf, 'PaperPosition', [0 0 6 5]);
    saveas(gcf, 'eqbcorrespondance_v10_eta', 'fig')
    saveas(gcf, 'eqbcorrespondance_v10_eta.eps', 'psc2')



