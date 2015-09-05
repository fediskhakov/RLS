%         Fedor Iskhakov, University Technology Sidney
%         John Rust, University of Maryland
%         Bertel Schjerning, University of Copenhagen
%         London Heathrow, August 2012

%% CLEAR MEMORY AND SWT SWITCHES
clc; clear; close all;
 
% Switches (Local to this script)
printparam=0;
allequilibria=1;
slideshow=0;
costshow=0;
simulate=0;
compile=1;
efficiency=1;
simulate=0;

%% COMPILE 
tic
if compile
    fprintf('COMPILING... '); 
    tic
%     mex -largeArrayDims leapfrog.c -DMAXEQB=3 -DPRINTeqbloop=1 -DPRINTeqbstr=0
    mex -largeArrayDims leapfrog.c -DMAXEQB=3 -DPRINTeqbstr=1
    tc=toc; fprintf('Compiled model in %1.10f (seconds)\n\n',tc);
end

%% SET BECNHMARK PARAMETERS
setup;

% Adjust mp structure - see setup for description
mp.eta=0;
mp.tpm=[0 1;
        1 0];

%mp.k1=8.3;
mp.k1=5;
%mp.k2=1;
mp.k2=0;
par.nC=5;                 % 4
capT=floor(1.5*par.nC);
mp.dt=5/par.nC;            % 5 is top cost
mp.df=exp(-0.05*mp.dt);
mp.c_tr=1;
mp.onestep=1;

% recalculate depend parameters
par.pti=f_pti(mp, par);
mp.nC=par.nC;  % add nC to mp

% Adjustments to sw structure - see setup for description
sw.alternate=true; %alternate move or simultanious move game
sw.esr=5; %equilibrium selection rule to be used: see setup.m or esr.c %ESR:always play first equilibrium
%sw.esr=1; %equilibrium selection rule to be used: see setup.m or esr.c %ESR:if available play mixed strategy
sw.esrmax=1000000; %run N feasible eqstrings (set equal to 192736405 for n=5)
sw.esrstart=0; %index of the first eqstring

% LOOK at particular equilibrium in simultanious move game
sw.esr=99;
sw.esrmax=1;
% sw.esrstart=940660880880; %monopoly MAXEQB=5
% sw.esrstart=195369009007; % effic 1.0009 MAXEQB=5
%sw.esrstart=516560655;
%sw.esrstart=47069207652; %mixed strategies everywhere in interior points
%sw.esrstart=47069204493;
%sw.esrstart=31467153051; %mixed strategy with efficient outcome
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
%mp_mon.c_tr=-1;
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

% mp_mon.tpm=flipud(fliplr(mp_mon.tpm)); %relable the firms for test of assymetry
tic; [bne_mon, br_mon, g_mon, eqbstr_mon]=leapfrog(par_mon,mp_mon,sw_mon); ts=toc;

fprintf('Solve model in %1.10f (seconds)\n\n\n\n',ts);

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

%Output efficiency measures in each state point
if efficiency
    for iC=1:numel(g)
        gmon=g_mon(iC).solution(~isnan(g_mon(iC).ec(:,1)),:);
        maxeqb=find(g(1).solution(:,2)==1,1,'first')-1;
        gmon=kron(gmon, ones(maxeqb,1));
        mask=~isnan(g(iC).ec(:,1)); %skip nans
        mask=mask & (g(iC).solution(:,4)<=g(iC).solution(:,5));
        %expected monopoly profits (=max social surplus) in every point
        monopoly_values=gmon(mask,11).*gmon(mask,8)+(1-gmon(mask,11)).*gmon(mask,7);
        if mp.nC<=15
            if sw.alternate
                fprintf('Level of the game %d\nCols: c1, c2, efficiency if m=1, if m=2\n',iC);
                disp([g(iC).solution(mask,[4 5]) g(iC).ec(mask,[9 10])./(monopoly_values*[1 1])]);
            else
                fprintf('Level of the game %d\nCols: c1, c2, ieqb, revenue-cost, moopoly profits, efficiency\n',iC);
                disp([g(iC).solution(mask,[4 5]) g(iC).solution(mask,3)+1 g(iC).ec(mask,9) monopoly_values g(iC).ec(mask,9)./monopoly_values]);
            end
        end
    end

    % eqbstr:
    % 1: lex index of
    % 2: number of repetions (sum of col2 numer of equilibria in total)
    % 3: v10
    % 4: v20 (x20 if alternating move)
    % 5: stat1 pure strategy dynamic equilibrium (0/1)
    % 6: stat2 symetric dynamic aquilibrium (0/1)
    % 7: stat3 leapfrogging eqb (high cost follower has poitsitive investment probability) (0/1)

    %% PRINT MONOPOLY AND DUOPOLY PROFITS
    fprintf('-------------------------------------------------------------------------------------\n');
    v10_mon=g_mon(par.nC).solution(1,7);
    v20_mon=g_mon(par.nC).solution(1,9);
    fprintf('Monopoly profits: %f\n', v10_mon+v20_mon);
    if sw.alternate
        fprintf('\nDuopoly  profits (alternating move)\n');
        v10_duo=g(par.nC).solution(1,7); 
        v20_duo=g(par.nC).solution(1,18);
    else
        fprintf('\nDuopoly  profits (simultaneous move)\n');
        v10_duo=g(par.nC).solution(1,7);
        v20_duo=g(par.nC).solution(1,9);
    end;
    fprintf(' Firm 1: %f\n', v10_duo);
    fprintf(' Firm 2: %f\n', v20_duo);
    fprintf(' Total : %f\n', v10_duo+v20_duo);
    fprintf('-------------------------------------------------------------------------------------\n');
    fprintf('\n\n');

    %eqbstr(:,8)=eqbstr(:,8)/v10_mon;
    %d=analyse(eqbstr);

end 


if simulate
    %% SIMULATE COST SEQUENCES
    % Simulate and plot realized sequences
    sp=simul.setup(par);         % Run setup for simulation module
    sp.T=capT;
    s=simul.sequence(g,sp, mp,par,sw.alternate);   % Sumulate sequences
    graph.CostSequence(s, mp, 'Duopoly equilibrium realization'); % plot sequences realized costs
    % graph.ProfitSequence(s, mp, 'Duopoly profits'); % plot sequences realized costs

    s_mon=simul.sequence(g_mon,sp, mp_mon,par_mon,1);   % Sumulate sequences
    graph.CostSequence(s_mon, mp, 'Monopoly equilibrium realization'); % plot sequences realized costs
    % graph.ProfitSequence(s_mon, mp_mon, 'Monpoly profits'); % plot sequences realized costs
    %return
end

%% All equilirbia
if allequilibria;
    if (mp.nC<=5 || (sw.alternate & mp.nC<=20));
        fprintf('SOLVE BENCHMARK MODEL - for all equilribria\n');
        sw.esr=99; %equilibrium selection rule to be used: see setup.m or esr.c
        sw.esrmax=1000000; %run N feasible eqstrings (set equal to 192736405 for n=5)
        sw.esrstart=0;  %

        tic;
        [bne, br, g, eqbstr, moncmp]=leapfrog(par,mp,sw,g_mon); %monopoly g structure for underinvestment checks
        [bne, br, g, eqbstr]=leapfrog(par,mp,sw);
        ts=toc;
        fprintf('Solved model in %1.10f (seconds)\n',ts);
        fprintf('doing maitnance to solution results\n');
        %eqbstr(:,1:min(10,size(eqbstr,2)));
        eqbstr(:,isnan(eqbstr(1,:)))=[];
        eqbstr=eqbstr';
        moncmp=moncmp';

        %Efficiency
        %eqbstr(:,8)=eqbstr(:,8)/(v10_mon+v20_mon);
        eqbstr(:,8)=eqbstr(:,8)/(v10_mon+v20_mon);

        if sw.alternate
            d=analyse(eqbstr,'alt');
            d=analyse(eqbstr,'notab','alt');
        else
            d=analyse(eqbstr);
            d=analyse(eqbstr,'notab');
        end

        %Graph and statistics
        desc='Pay-off map for all found equilibria';
        [a, b]=graph.EqbstrPlot(eqbstr,[2 8],[v10_mon+v20_mon],mp, sw,desc);

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
        fprintf('      Lex index   # of ESR        v10    v20/x20       pure     syemmetric  leapfrog   efficiency\n');
        fprintf('---------------------------------------------------------------------------------------------------\n');
        fprintf('%15d %10d %10g %10g %10d %10d %10d %15g \n', eqbstr_mon_f1(:,1:8)')
        fprintf('%15d %10d %10g %10g %10d %10d %10d %15g \n', eqbstr_mon_f2(:,1:8)')
        fprintf('\n');
        fprintf('Number of monopoly equilirbia   : %f\n', sum(eqbstr_mon_f1(:,2)+eqbstr_mon_f2(:,2)));
        fprintf('Number of equilirbia            : %f\n', sum(eqbstr(:,2)));
        fprintf('Monopoly equilibria (share of all equilirbia) : %f\n', sum(eqbstr_mon_f1(:,2)+eqbstr_mon_f2(:,2))/sum(eqbstr(:,2)));
        fprintf('---------------------------------------------------------------------------------------------------\n\n');

        %MORE Graphs
        %for the figures in the paper:
        %skip long title by passing [] in place of mp
        graph.EqbstrPlot(eqbstr,[5 8],(v10_mon+v20_mon),mp,sw,desc); %pure
        %graph.EqbstrPlot(eqbstr,[6 8],(v10_mon+v20_mon),[],sw,desc); %symmetric
        %graph.EqbstrPlot(eqbstr,[7 8],(v10_mon+v20_mon),[],sw, desc);%leapfrog

        d=analyse(eqbstr);

        dummystr = '[((eqbstr(:,3)~=0).*(eqbstr(:,4)~=0)) ((eqbstr(:,3)==0) | (eqbstr(:,4)==0))]';
        graphstring='graph.EqbstrCDFPlot(eqbstr,[0], dummystr, lbl ,mp,sw,''Distribution of found equilibria'',0.5, mfig)';
        graph.EqbstrCDFPlot(eqbstr,[0 -5 6 7], [], [] ,[],sw,'Distribution of found equilibria',1);
        %movie
        % dr ='../../../movies/';
        % cdfmovie(par,mp,sw,'mp.k1',[0 20],[dr 'CDF.k1.n4_tmp'], .5, dummystr, {'vi ne 0','vi=0'} , graphstring);
      
    end;
end

return


%%% OLD STUFF BELOW

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
