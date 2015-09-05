%  Results.m: generates graphs for section 4 in Iskhakov, Rust and Schjerning (2012)
%         Fedor Iskhakov, University of New South Wales
%         John Rust, University of Maryland
%         Bertel Schjerning, University of Copenhagen

%% CLEAR MEMORY AND SWT SWITCHES
clc; clear; close all;

% to be modified by user 
paperdir='/Users/bertelschjerning/Dropbox/IRS/leapfrogging/_IOpaper/current/0pdf/';

%%Switches
compile=1;

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
mp.c_tr=-1; 
mp.onestep=1;
mp.tpm=[0 1; 1 0];
mp.k2=0; 
par.nC=25;                 
mp.nC=par.nC;  % add nC to mp







%%%%%%%%%%%%%%%%%%%
%% FIGURE 4
%%%%%%%%%%%%%%%%%%%
mp.k1=5;
for i=1:2;
	mpvec(i)=mp; 
	swvec(i)=sw; 
end

i=1; mpvec(i).c_tr=1;  										
i=2; mpvec(i)=mp; swvec(i).esr=7; swvec(i).alternate=false;	

tit1(1)={'Panel (a): Non-monotonic tech. progress'}
tit1(2)={'Panel (b): Simultaneous move'}


tit2(1)={'Panel (c): Non-monotonic tech. progress'}
tit2(2)={'Panel (d): Simultaneous move'}


for i=1:numel(mpvec);
	par.nC=5;                 
    if swvec(i).alternate==false		
    	par.nC=5;                 
	end
	mpvec(i).nC=par.nC;

	% recalculate depend parameters
	mpvec(i).dt=5/par.nC;            % 5 is top cost
	mpvec(i).df=exp(-0.05*mpvec(i).dt);
	par.pti=f_pti(mpvec(i), par);

  	swvec(i).esr=99; %equilibrium selection rule to be used: see setup.m or esr.c
    swvec(i).esrmax=10000; %run N feasible eqstrings (set equal to 192736405 for n=5)
    swvec(i).esrstart=0;  %

	%% SET MONOPOLY PARAMETERS
	mp_mon=mpvec(i);
	par_mon=par;
	sw_mon=swvec(i);
	mp_mon.tpm=[1 0; 1 0];
	sw_mon.alternate=true; % always  solve alternate move to obtain monopoly profits
	sw_mon.esr=1; %equil5ibrium selection rule to be used: see setup.m or esr.c
	sw_mon.esrmax=1; %run N feasible eqstrings (set equal to 192736405 for n=5)

	%% SOLVE MONOPOLY MODEL
	[bne_mon, br_mon, g_mon, eqbstr_mon]=leapfrog(par_mon,mp_mon,sw_mon); 

	v10_mon=g_mon(par.nC).solution(1,7);
	v20_mon=g_mon(par.nC).solution(1,9);

	%% SOLVE MODEL AT PARAMETRS mpvec(i), par, swvec(i)
	[bne, br, g, eqbstr, moncmp]=leapfrog(par,mpvec(i),swvec(i),g_mon);

	eqbstr(:,isnan(eqbstr(1,:)))=[];
	eqbstr=eqbstr';
	moncmp=moncmp';

    %Efficiency
    %eqbstr(:,8)=eqbstr(:,8)/(v10_mon+v20_mon);
    eqbstr(:,8)=eqbstr(:,8)/(v10_mon+v20_mon);

    [a, b]=graph.EqbstrPlot(eqbstr,[2 8],(v10_mon+v20_mon),[],sw,tit1{i}); 
%    title(tit1{i});

    set(gcf, 'PaperPosition', [0 0 6 6]);
	saveas(gcf, sprintf('%sTri_%d',paperdir,  i), 'fig')
	saveas(gcf, sprintf('%sTri_%d.eps',paperdir,  i), 'psc2')


%   d=analyse(eqbstr);

    dummystr = '[((eqbstr(:,3)~=0).*(eqbstr(:,4)~=0)) ((eqbstr(:,3)==0) | (eqbstr(:,4)==0))]';
    graphstring='graph.EqbstrCDFPlot(eqbstr,[0], dummystr, lbl ,mp,sw,tit2{i},0.5, mfig)';
  %  graph.EqbstrCDFPlot(eqbstr,[0 -5 6 7], dummystr, {'vi ne 0','vi=0'} ,[],sw,tit2{i},1);
    graph.EqbstrCDFPlot(eqbstr,[0 7 -5 ], [], [] ,[],sw,tit2{i},1);
 %   title(tit2{i});
    set(gcf, 'PaperPosition', [0 0 6 6]);
	saveas(gcf, sprintf('%sCDF_%d',paperdir,  i), 'fig')
	saveas(gcf, sprintf('%sCDF_%d.eps',paperdir,  i), 'psc2')


end





%%%%%%%%%%%%%%%%%%%
%% FIGURE 6 and 7
%%%%%%%%%%%%%%%%%%%
clear mpvec, swvec;
mp.k1=2;  
for i=1:7;
	mpvec(i)=mp; 
	swvec(i)=sw; 
end
	
% Adjust mp structure - see setup for description
% FIGURE 6
i=1; mpvec(i).k1=2;  mpvec(i).nC=100; swvec(i).esr=99;		    tit(i)={'Panel (a): Preemption and rent-dissipation'}
i=2; mpvec(i).k1=2;  mpvec(i).nC=25;  swvec(i).esr=99;		    tit(i)={'Panel (b): Underinvestment'}
i=3; mpvec(i).k1=.5; mpvec(i).nC=25;  swvec(i).esr=99;		    tit(i)={'Panel (c): Leap-frogging'}


% FIGURE 7
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
	sw_mon.esrmax=1; %run N feasible eqstrings (set equal to 192736405 for n=5)

	%% SOLVE MONOPOLY MODEL
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



