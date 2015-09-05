%% dt-->0 in deterministic technological progress deterministic alternating move

%% VERY SLOW

%         Fedor Iskhakov, University Technology Sidney
%         John Rust, University of Maryland
%         Bertel Schjerning, University of Copenhagen
%         Sidney, March 2012

close all;

%% COMPILE 
if true
    fprintf('Compiling... '); 
    tic
    mex -largeArrayDims leapfrog.c -DMAXEQB=3 -DPRINTeqbloop=0 -DPRINTeqbstr=0
    tc=toc; fprintf('Compiled model in %1.10f (seconds)\n\n',tc);
end

%% PARAMETERS
setup;
mp.dt=1; % initial value
mp.k1=10.0;
mp.k2=1.0;
mp.eta=0;
mp.c_tr=-0.5;%deterministic technological progress!!!

mp.onestep=1;  
mp.tpm=[0 1;
        1 0];%deterministic alternating move
sw.alternate=true; %alternate move
sw.esr=99; %first equilibrium every time - we know there is only 1 everywhere


%% LOOP OVER dt
dtdata=[];
for step=1:10
    %dt (bisections)
    mp.dt=.5^step;

    % recalculate depend parameters
    [par mp]=f_update_params(par,mp);

    fprintf('dt-->0 step %d: dt=%1.3e nC=%3d : ',step,mp.dt,mp.nC);

    %MONOPOLY SOLUTION
    mpm=mp;
    mpm.tpm=[1 0; 1 0];
    swm.alternate=true;
    swm.analytical=false;
    swm.esr=5;
    swm.esrmax=10;
    swm.esrstart=0;
    fprintf('monopoly ');
    tic;
    [a,b,gm,c]=leapfrog(par,mpm,swm);
    tsm=toc;
    fprintf('(%1.10fs) : ',tsm);
    mon=gm(par.nC).solution(1,7)+gm(par.nC).solution(1,9);

    fprintf('solver ');
    tic; 
    [bne, br, g, eqbstr]=leapfrog(par,mp,sw);
    ts=toc;
    eqbstr(:,isnan(eqbstr(1,:)))=[];
    eqbstr=eqbstr';
    fprintf('(%1.10fs) : ',ts);

    fprintf('nrEqb=%d df=%1.5g Mon=%g v10=%g x10=%g \n',sum(eqbstr(:,2)),mp.df,mon,g(par.nC).solution(1,7),g(par.nC).solution(1,16));

    %save results
    dtdatalabels={'dt' 'value_firm1' 'value_monopoly' 'n' 'df' 'number of equilibria' 'runtime_firm1' 'runtime_monopoly'};
    dtdata=[dtdata; mp.dt max(g(par.nC).solution(1,16),g(par.nC).solution(1,7)) mon mp.nC mp.df sum(eqbstr(:,2)) ts tsm]; %use the higher value between my turn and not my turn

    %save to dist
    save('dtdata','dtdata','dtdatalabels');
    
end


