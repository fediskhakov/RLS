%         Fedor Iskhakov, University Technology Sidney
%         John Rust, University of Maryland
%         Bertel Schjerning, University of Copenhagen
%         Sidney, March 2012

%% CLEAR MEMORY AND SWT SWITCHES
clc
clear;
close all;
 
% Switches (Local to this script)
allequilibria=0;
slideshow=0;
costshow=0;
simulate=1;
compile=1;
efficiency=0;

%% COMPILE 
tic
if compile
    fprintf('COMPILING... '); 
    tic
%     mex -largeArrayDims leapfrog.c -DMAXEQB=3 -DPRINTeqbloop=1 -DPRINTeqbstr=0
    mex -largeArrayDims leapfrog.c -DMAXEQB=5
    tc=toc; fprintf('Compiled model in %1.10f (seconds)\n\n',tc);
end

%% SET BECNHMARK PARAMETERS 
setup;
% Adjust mp structure - see setup for description
% mp.k1=10; mp.k2=1.0; mp.c_tr=1;  % ASSYMMETRIC
% mp.k1=1; mp.k2=10; mp.c_tr=0.05; % SYMMETRIC
mp.eta=0;
mp.tpm=[0 1;
        1 0]; 
% par.nC=100;
% mp.dt=80/par.nC;
% mp.df=mp.df/par.nC;
    
mp.c_tr=-.5;
mp.onestep=1;  
mp.beta_a=1.8; 
mp.beta_b=0.4;    


% mp.tpm=[0 1;
%         1 0]; 

% mp.tpm=[0.5 0.5;
%         0.5 0.5]; 

% mp.k1=40; mp.k2=1; mp.c_tr=-0.5; mp.eta=0; 
% mp.k1=10; mp.k2=1; mp.c_tr=0.05; mp.eta=.0; 
% mp.k1=10; mp.k2=1.0; mp.c_tr=1;  % ASSYMMETRIC

% Adjustments to par structure - see setup for description
par.nC=4/mp.dt;
par.pti=f_pti(mp, par);

% recalculate depend parameters
[par mp]=f_update_params(par,mp);

% Adjustments to sw structure - see setup for description
% sw.analytical=false; 
sw.alternate=false; %alternate move or simultanious move game
sw.esr=99; %equilibrium selection rule to be used: see setup.m or esr.c
sw.esrmax=1000000; %run N feasible eqstrings (set equal to 192736405 for n=5)
sw.esrstart=0; %index of the first eqstring

%% SET MONOPOLY PARAMETERS 
mp_mon=mp;
par_mon=par;
sw_mon=sw;

% Adjustments to par structure - see setup for description
mp_mon.tpm=[1 0; 1 0];

% Adjustments to sw structure - see setup for description
sw_mon.alternate=true; % always  solve alternate move to obtain monopoly profits
sw_mon.esr=1; %equilibrium selection rule to be used: see setup.m or esr.c
sw_mon.esrmax=1; %run N feasible eqstrings (set equal to 192736405 for n=5)

%% SOLVE MONOPOLY MODEL
fprintf('SOLVE MONOPOLY MODEL\n');
% set parameters to benchmark parameters
fprintf('Parameters are:\n');
mp_mon
tpm=mp_mon.tpm
par_mon
sw_mon

fprintf('Solving model... \n')
tic; [bne_mon, br_mon, g_mon, eqbstr_mon]=leapfrog(par_mon,mp_mon,sw_mon); ts=toc;
fprintf('Solve model in %1.10f (seconds)\n\n',ts);

v10=g_mon(par.nC).solution(1,7);
v11=g_mon(par.nC).solution(1,8);
v20=g_mon(par.nC).solution(1,9);
v21=g_mon(par.nC).solution(1,10);

mon=g_mon(par.nC).solution(1,7)+g_mon(par.nC).solution(1,9);


%% SOLVE MODEL AT PARAMETRS mp, par, sw
if allequilibria
    fprintf('SOLVE BENCHMARK MODEL\n');
    fprintf('Parameters are:\n');
    mp
    tpm=mp.tpm
    par
    sw.analytical=0
    
    fprintf('Solving model... \n')
    tic;
    [bne, br, g, eqbstr]=leapfrog(par,mp,sw);
    ts=toc;
    fprintf('Solved model in %1.10f (seconds)\n',ts);
    fprintf('doing maitnance to solution results\n');
    eqbstr(:,1:10)
    eqbstr(:,isnan(eqbstr(1,:)))=[];
    eqbstr=eqbstr';
end

%% GRAPHICS
if allequilibria
    [a, b]=graph.EqbstrPlot(eqbstr,[2 5],[v10+v20],mp, sw,'All eqb')
end% eqbstr:
% 1: lex index of
% 2: number of repetions (sum of col2 numer of equilibria in total)
% 3: v10
% 4: v20 (x20 if alternating move)
% 5: stat1 pure strategy dynamic equilibrium (0/1)
% 6: stat2 symetric dynamic aquilibrium (0/1)
% 7: stat3 leapfrogging eqb (high cost follower has poitsitive investment probability) (0/1)

%% EFFICIENCY    fprintf('Parameters:\n')
if efficiency

    disp(par);
    disp(mp);
    disp(sw);
    fprintf('Transition probability for m:\n')
    disp(mp.tpm);
    
    tic; 
    [bne, br, g, eqbstr]=leapfrog(par,mp,sw);
    ts=toc;
    fprintf('Solved model in %1.10f (seconds)\n',ts);
    eqbstr(:,isnan(eqbstr(1,:)))=[];
    eqbstr=eqbstr';
    return;

    %Efficiency
    %eqbstr(:,8)=1-eqbstr(:,8)/mon;
    eqbstr(:,8)=eqbstr(:,8)/mon;

    %eqbsum=graph.EqbstrPlot (eqbstr,[2 2],mon,mp,sw,'Pay-off map for all found equilibria');
    %fprintf('%35s : %5d\n','Total number of equilibria',eqbsum(1), ...
    %                       'Number of distinct equilibria',eqbsum(2), ...
    %                       'Number of pure strategy equilibria',eqbsum(3), ...
    %                       'Number of symmetric equilibria',eqbsum(4), ...
    %                       'Number of leapfrogging equilibria',eqbsum(5));

desc='';
eqbstr(:,8)
%for the figures in the paper:
    %skip long title by passing [] in place of mp
    graph.EqbstrPlot(eqbstr,[8 8],mon,mp,sw,desc);
    %graph.EqbstrPlot(eqbstr,2,mon,[],sw,desc);
    %graph.EqbstrPlot(eqbstr,[5 2],mon,[],sw,['(2) ' desc ': pure equilibria']);
    %graph.EqbstrPlot(eqbstr,[6 2],mon,[],sw,['(3) ' desc ': symmetric equilibria']);

end

%% DISPLAY SOLUTION FOR A SINGLE EQUILIBRIUM SELECTION RULE
% CHANGE BENCHMARK PARAMETERS HERE
sw.esr=1; % ESR_HighCostInvestsMix   =3;
jC=1;
% Adjustments to par structure - see setup for description
mp.dt=1;
mp.df=exp(-0.05*mp.dt);

par.nC=4/mp.dt;
par.pti=f_pti(mp, par);

scrsz = get(0,'ScreenSize');
figure1 = figure('Name','Endgame Solution','Color',[1 1 1],'Position',[1 1 scrsz(3) scrsz(4)]);
[bne, br, g, eqbstr]=leapfrog(par,mp,sw);
graph.solution(g,jC, par,1,figure1);

gcost=g;
gcost(jC).solution(:,7:10)=g(jC).ec(:,1:4);
figure2 = figure('Name','Endgame Expected Cost','Color',[1 1 1],'Position',[1 1 scrsz(3) scrsz(4)]);
graph.solution(gcost,jC, par,1,figure2);

if  sw.alternate==true
    figure2 = figure('Name','Endgame Solution by turn','Color',[1 1 1],'Position',[1 1 scrsz(3) scrsz(4)]);
    graph.solution(g,1, par,2,figure2);
    % NOTE: when second to last argument (plottype==2) firm2 value functions 
    % and investment probabilities  in the lower panel are repalced by 
    % values and probabilites of firm 2 when it is not it's turn to invest
    % (only valid for alternating move game)
    % Last argument is a optional figure handle
end

%% DISPLAY SOLUTION SLIDESHOW FOR A SINGLE EQUILIBRIUM (plottype=1)
if slideshow
    scrsz = get(0,'ScreenSize');
    figure3 = figure('Name','Endgame Solution by turn','Color',[1 1 1],'Position',[1 1 scrsz(3) scrsz(4)]);
    for iC=1:par.nC-1
        pause(10/par.nC);
        graph.solution(g,iC, par,1,figure3);
    end
end
symmetry_checker(g)

%% DISPLAY SOLUTION SLIDESHOW FOR A SINGLE EQUILIBRIUM (plottype=1)
if costshow
    figure3 = figure('Name','Expected costs','Color',[1 1 1],'Position',[1 1 scrsz(3) scrsz(4)]);
    gcost=g;
    for iC=2:par.nC-1
        gcost(iC).solution(:,7:10)=g(iC).ec(:,1:4);
%         pause(5/par.nC);
        graph.solution(gcost,iC, par,1,figure3);
    end
end
%% DISPLAY SOLUTION SLIDESHOW FOR A SINGLE EQUILIBRIUM (plottype=2)
    % NOTE here plottype==2 such that firm2 value functions 
    % and investment probabilities  in the lower panel are repalced by 
    % values and probabilites of firm 2 when it is not it's turn to invest
    % (only to run for alternating move game)
    
if slideshow && sw.alternate==true
    figure3 = figure('Name','Model solution - Slideshow','Color',[1 1 1],'Position',[1 1 scrsz(3) scrsz(4)]);
    pause;
    for iC=1:par.nC-1;
        pause(0.1);
        graph.solution(g,iC, par,2, figure3);
    end
end
% return


%% SECTION 5.0: Simulate and plot realized sequences
if simulate;
    sp=simul.setup(par)          % Run setup for simulation module
    s=simul.sequence(g,sp, mp,par,1)   % Sumulate sequences
    graph.CostSequence(s, mp, 'Equilibrium realization'); % plot sequences realized costs
    graph.ProfitSequence(s, mp, 'Profits'); % plot sequences realized costs
end

return






























%% --------------------------------------------------------------------------------------------
%  ACHTUNG! CODE BELOW THIS LINE IS NOT YET FUNCTIONAL
% --------------------------------------------------------------------------------------------


%% SECTION 4.1: BEST RESPONSE FUNCTIONS
iC=1;
spacing=0.001;
figure('Name','Best response function');
pvec=0.0000001:.001:.998;
[bne, br ,g]=leapfrog(par,mp);
graph.br(br, g, mp.eta, c1,c2,iC, spacing);
%% SECTION 4.2: SECOND ORDER BEST RESPONSE FUNCTION COLARED BY ROOTS
iF=1;                 % Plot for firm If
figure('Name','Second order best response function - FIRM 1');
graph.br2byroots(br, g, mp.eta, c1,c2,iC,iF, 0.001);

%% SECTION 4.3: BEST RESPONSE FUNCTION SLIDE SHOW
iC=1;
if slideshow;
eta0=mp.eta;
    figure('Name','Best response function');
    for eta=0:0.1:1;
        mp.eta=eta;
        [bne, br ,g]=leapfrog(par,mp);
        hold on;
        graph.br(br, g, eta, c1,c2,iC, spacing);
        pause(.05);
        hold off
    end;
mp.eta=eta0;
[bne, br ,g]=leapfrog(par,mp);
end
%% SECTION 4.4: SECOND ORDER BEST RESPONSE FUNCTION SLIDE SHOW
if slideshow;
eta0=mp.eta;
    iC=1; iF=1;
    figure('Name','Second order best response function - Slideshow');
    for eta=0:0.1:1;
        mp.eta=eta;
        [bne, br ,g]=leapfrog(par,mp);
        % hold on;
        graph.br2(br, g, eta, c1,c2,iC,iF, spacing);
        pause(.2);
        hold off
    end
mp.eta=eta0;
[bne, br ,g]=leapfrog(par,mp);
end
%% TESTING AREA:
           
%% fixed point algorithm applied on test function
% fplot(@fun.testfun,[0,1])
figure('Name','Fixed point algorithm applied on test function');
lo=0; hi=1;
y=lo:.0001:hi;
x=fun.testfun(y);
hold on;
plot(x,y);
hold on
plot(lo:0.001:hi,lo:0.001:hi, '-r')
    
fx= @(x)fun.testfun(x);
par.maxit=2000;
par.ctol=1e-10;
eqbinfo=fun.findeqb(par, lo, hi, fx);
for i=1:size(eqbinfo,1)
    line([eqbinfo(i,4); eqbinfo(i,4)], [lo; hi],'Color', 'r', 'LineStyle', ':');
    line([eqbinfo(i,2); eqbinfo(i,2)], [lo; hi], 'Color', 'b','LineStyle', ':');
end


%% fixed point algorithm applied on test function
% plot inverse of inverse
lo=0; hi=1;
clf;
hold off;
y=lo:.0001:hi;
x=fun.testfun(y);
plot(y,x);
hold on
plot(lo:0.001:hi,lo:0.001:hi, '-r')

fx= @(x)fun.testfun(x);
par.maxit=2000;
par.ctol=1e-10;
eqbinfo=fun.findeqb(par, lo, hi, fx);
for i=1:size(eqbinfo,1)
    line([lo; hi],[eqbinfo(i,4); eqbinfo(i,4)],'Color', 'r', 'LineStyle', ':');
    line([lo; hi], [eqbinfo(i,2); eqbinfo(i,2)], 'Color', 'b','LineStyle', ':');
end

%% Check of findeqb and FindMinMaxBr2
% select polynomial coefficient for combination of cost's chosen in section 1.2 above
brc=br(iC).br;
i= brc(:,3) == c1 & brc(:,4) == c2;
brc=brc(i,:);

% Find eqb for a given root combination as fixed points on second order best response function 
r=4;% root combination to look at (takes values 1-4)
[lo, hi,invbr2_min, invbr2_max]=fun.FindMinMaxBr2(brc, mp,par, 0.0000001, .9999999, iF, r);
fx= @(x)fun.invbr2byroot(brc, mp.eta, x, iF, r);
eqbinfo=fun.findeqb(par, lo, hi, fx);

% CHECK OF FindEqb 
hold on
for i=1:size(eqbinfo,1) 
    line([eqbinfo(i,4); eqbinfo(i,4)],[0; 1],'Color', 'r', 'LineStyle', ':');
%    line([eqbinfo(i,2); eqbinfo(i,2)],[0; 1], 'Color', 'b','LineStyle', ':');
end
for i=2:size(eqbinfo,1)-1
    line([eqbinfo(i,2); eqbinfo(i,2)],[0; 1], 'Color', 'b','LineStyle', ':');
end
return

% CHECK OF FindMinMaxBr2
brc=br(iC).br;
c=g(iC).c;
i=find(brc(:,3) == c1 & brc(:,4) == c2);
brc=brc(i,:);
pimax=0.999999999;
pimin=0.000000001;

for r=1:4;
    [pimin_out(r), pimax_out(r),invbr2_min, invbr2_max]=fun.FindMinMaxBr2(brc, mp,par, pimin, pimax, iF, r)
     eqb=fun.FindEqb(br, par, mp, c1,c2,iC,1);
end
return
