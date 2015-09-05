%         Fedor Iskhakov, University Technology Sidney
%         John Rust, University of Maryland
%         Bertel Schjerning, University of Copenhagen
%         November 2011

%% SECTION 1.1: CLEAR MEMORY
clc
clear;
close all


%% SECTION 1.2: SET PARAMETERS
setup; 
% Adjustments to par structure - see setup for description
% par.ctol=0.0000001;
par.nC=11;  
par.nP=200;

% Adjust mp structure - see setup for description
mp.eta=0; 
mp.sigma=0;

slideshow=0;
simulate=1;
%% SECTION 2.1: SOLVE MODEL IN SERIAL 
% mex -largeArrayDims leapfrog.c 
nrep=1;    %% increase to increase time mesurement precision when cpu time is tiny
tic
for i=1:nrep;
    [bne, br ,g]=leapfrog(par,mp);
end
ts=toc/nrep;
fprintf('Running time to solve model (serial code): %1.10f seconds\n',ts);
a=g(2).solution;
% save(['testers' filesep 'JRvsC' filesep 'g.mat'],'g');

%% SECTION 2.1: SOLVE MODEL IN PARALLEL (XXX NOT ACTIVE)
% mex -largeArrayDims COMPFLAGS="$COMPFLAGS -openmp" LDFLAGS="$LDFLAGS -openmp"  leapfrog.c 
% tic
% for i=1:1;
% [bne, br ,g]=leapfrog(par,mp);
% end;
% tp=toc/nrep;

% fprintf('Running time to solve model (parallel code): %1.10f seconds\n',tp);
% fprintf('Running time to solve model (serial code): %1.10f seconds\n',ts);
% fprintf('Speed factor (serial/parallel): %1.10f seconds\n',ts/tp);


%% SECTION 3.1: DISPLAY SOLUTION
figure('Name','Model solution');
graph.solution(g,1, par);
figure('Name','Model solution - Edges');
graph.edges(g,1, par);

%% SECTION 2.1: SLIDESHOW OF SOLUTION
if slideshow
figure('Name','Model solution - Slideshow');
    pause(1);
    for iC=1:par.nC-1;
        graph.solution(g,iC, par);
        pause(.1);
    end
end

%% SECTION 5.0: Simulate and plot realized sequences
if simulate;
    sp=simul.setup(par)          % Run setup for simulation module
    s=simul.sequence(g,sp, mp)   % Sumulate sequences
    graph.CostSequence(s, mp, 'Equilibrium realization'); % plot sequences realized costs
    graph.ProfitSequence(s, mp, 'Profits'); % plot sequences realized costs
end
%% SECTION 4.0.1: SET PARAMETERS FOR BEST RESPONSE FUNCTIONS
close all
clc
spacing=0.001;
mp.eta=0; c1=1; c2=1; iC=1; 
mp.eta=.3; c1=1.5; c2=0.5; iC=1; 
eta=0.3; c1=1; c2=1; iC=1; 
eta=0.3; c1=3.5; c2=0.5; iC=1; 
eta=0.3; c1=4.5; c2=1.5; iC=1; 
eta=0.3; c1=1.5; c2=0.5; iC=1; 
eta=0.5; c1=1.5; c2=1.5; iC=1;  %% one unstabel fixed points
eta=0.3; c1=1; c2=1; iC=1;  %% three fixed poinst
eta=0.3; c1=1.5; c2=4.5; iC=1;  %% two fixed points
eta=0.3; c1=2.5; c2=2.5; iC=1;  %% one unstabel fixed points
eta=0.15; c1=3.5; c2=3; iC=5; 
par.ctol=1e-12;
mp.eta=eta;
iF=1;

% SECTION 4.0.2: BEST RESPONSE FUNCTION
figure('Name','Best response function');
pvec=0.0000001:.001:.998;
[bne, br ,g]=leapfrog(par,mp);
graph.br(br, g, mp.eta, c1,c2,iC, spacing);

% % SECTION 6.0: SECOND ORDER BEST RESPONSE FUNCTION
% figure('Name','Second order best response function - FIRM 1');
% graph.br2(br, g, mp.eta, c1,c2,iC,1, spacing);
% 
% figure('Name','Second order best response function - FIRM 2');
% graph.br2(br, g, mp.eta, c1,c2,iC,2, spacing);

% SECTION 6.0: SECOND ORDER BEST RESPONSE FUNCTION COLARED BY ROOTS

figure('Name','Second order best response function - FIRM 1');
    [bne, br ,g]=leapfrog(par,mp);
    graph.br2byroots(br, g, eta, c1,c2,iC,iF, 0.001);
    pause(0.1);

brc=br(iC).br;
i= brc(:,3) == c1 & brc(:,4) == c2;
brc=brc(i,:);

% Find eqb as fixed points on second order best response function 
r=4;
[lo, hi,invbr2_min, invbr2_max]=fun.FindMinMaxBr2(brc, mp,par, 0.0000001, .9999999, iF, r);
fx= @(x)fun.invbr2byroot(brc, eta, x, iF, r);
clc;
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
r=2;
brc=br(iC).br;
c=g(iC).c;
i=find(brc(:,3) == c1 & brc(:,4) == c2);
brc=brc(i,:);
pimax=0.999999999;
pimin=0.000000001;

for r=1:4;
    [pimin_out(r), pimax_out(r),invbr2_min, invbr2_max]=fun.FindMinMaxBr2(brc, mp,par, pimin, pimax, iF, r)
end
return

eqb=fun.FindEqb(br, par, mp, c1,c2,iC,1);


%%
for r=1:4;
[invbr2,domain_id]=fun.invbr2byroot(brc, eta, p_i, iF, r);
end

%% COMPARE FindEqb with C-program (works only for eta=1)
i=find(g(iC).solution(:,4) == c1 & g(iC).solution(:,5) == c2);
load g_eta1.mat; % load solution for eta=1
costvector=g(iC).solution(i,4:5)
eqbvector=g(iC).solution(i,11:12)

%% test invbr2byroot(brc, eta, p_i, iF, -1, -1);
brc=br(iC).br;
c=g(iC).c;
i=find(brc(:,3) == c1 & brc(:,4) == c2);
brc=brc(i,:);
p_i=0.5;
iF=1;
clc
for ir=1:4;
    [root(ir), id(ir)]=fun.invbr2byroot(brc, eta, p_i, iF, ir);
end
[root; id]
allroots=fun.invbr2(brc, eta, p_i, iF)


%% SECTION 5.0: BEST RESPONSE FUNCTION SLIDE SHOW
iC=1;
if slideshow;
    figure('Name','Best response function');
    for eta=0:0.1:1;
        mp.eta=eta;
        [bne, br ,g]=leapfrog(par,mp);
        hold on;
        graph.br(br, g, eta, c1,c2,iC, spacing);
        pause(.2);
        hold off
    end;
end
%% SECTION 5.0: SECOND ORDER BEST RESPONSE FUNCTION SLIDE SHOW
if slideshow;
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
end

           
%% test function
% fplot(@fun.testfun,[0,1])
close all
lo=0; hi=1;
y=lo:.0001:hi;
x=fun.testfun(y);
hold on;
plot(x,y);
hold on
plot(lo:0.001:hi,lo:0.001:hi, '-r')
    
fx= @(x)fun.testfun(x);
clc;
par.maxit=2000;
par.ctol=1e-14;
eqbinfo=fun.findeqb(par, lo, hi, fx);
for i=1:size(eqbinfo,1)
    line([eqbinfo(i,4); eqbinfo(i,4)], [lo; hi],'Color', 'r', 'LineStyle', ':');
    line([eqbinfo(i,2); eqbinfo(i,2)], [lo; hi], 'Color', 'b','LineStyle', ':');
end

%% // plot inverse of inverse
lo=0; hi=1;
clf;
hold off;
y=lo:.0001:hi;
x=fun.testfun(y);
plot(y,x);
hold on
plot(lo:0.001:hi,lo:0.001:hi, '-r')

fx= @(x)fun.testfun(x);
clc;
par.maxit=2000;
par.ctol=1e-10;
eqbinfo=fun.findeqb(par, lo, hi, fx);
for i=1:size(eqbinfo,1)
    line([lo; hi],[eqbinfo(i,4); eqbinfo(i,4)],'Color', 'r', 'LineStyle', ':');
    line([lo; hi], [eqbinfo(i,2); eqbinfo(i,2)], 'Color', 'b','LineStyle', ':');
end



