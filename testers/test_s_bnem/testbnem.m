%         Fedor Iskhakov, University Technology Sidney
%         John Rust, University of Maryland
%         Bertel Schjerning, University of Copenhagen
%         March 2011


% to be dealt with: 
% bne.m does not appropriately take outside good into account when computing porfits
% bnem.c does not appropriately take outside good into account when sigam>0

%% SECTION 1: CLEAR MEMORY
clear
clc
close all

%% SECTION 2: SET PARAMETRS
setup; 
% Adjustments to par structure - see setup for description
par.ctol=0.000001;
par.nC=10;  
par.nP=60;

% Adjust mp structure - see setup for description
% mp.k1=2; 
mp.eta=0; 
mp.sigma=0.0000001;
mp.c_og=3;
par.og=0;    

%% SECTION 3a: SOLVE B_N equilibrium using bnem.c compute in parallel
mex -largeArrayDims COMPFLAGS="$COMPFLAGS -openmp" LDFLAGS="$LDFLAGS -openmp" bnem.c

tic
vbnem=bnem(par,mp);
tap=toc;
fprintf('Running time to solve model in C in Parallel: %1.10f seconds\n',tap);

%% SECTION 3b: SOLVE B_N equilibrium using bnem.c  - compute serially
mex -largeArrayDims bnem.c

tic
vbnem=bnem(par,mp);
ta=toc;
fprintf('Running time to solve model in C - serial compilation: %1.10f seconds\n',ta);
fprintf('Speed factor in parallel: %1.10f seconds\n',ta/tap);

%% SECTION 3: SOLVE B_N equilibrium using bne.m
clc;
tic
cgrid=NaN(par.nC,1);
for i=0:(par.nC-1);
    cgrid(i+1)=par.cmin+i*(par.cmax-par.cmin)/(par.nC-1);
end

vbne=NaN(par.nC*par.nC,10);
for ic1=0:par.nC-1;
for ic2=0:par.nC-1;
    i=ic1*par.nC+ic2;
vbne(i+1,1)=ic1;
vbne(i+1,2)=ic2;
vbne(i+1,3)=cgrid(ic1+1);
vbne(i+1,4)=cgrid(ic2+1);

bneinfo=bne(cgrid(ic1+1),cgrid(ic2+1),par.og, mp.c_og,mp.sigma);
vbne(i+1,5:10)=bneinfo;     % profits, firm 1
end
end
tb=toc;

%% SECTION 3: COMPARE

checkmat1=vbne-vbnem;
% fprintf('        c1        c2   p(1)_JR   p(2)_JR    p(1)_c    p(2)_c  pf(1)_JR  pf(2)_JR   pf(1)_c   pf(2)_c\n');
% format short
% display([vbne(:,3:4) vbne(:,5:6) vbnem(:,5:6)   vbne(:,7:8) vbnem(:,7:8)]);

fprintf('MAX(ABS(DIFF))\n', max(abs(checkmat1)));
fprintf('     p(1)      p(2)      pf(1)      pf(2)       s(1)       s(2)\n %f   %f   %f   %f   %f   %f%', max(abs(checkmat1)));
fprintf('\n');



%% SECTION 3: COMPARE CPU TIME
fprintf('Running time to solve model in C in Parallel: %1.10f seconds\n',tap);
fprintf('Running time to solve model in C - serial compilation: %1.10f seconds\n',ta);
fprintf('Running time to solve model in MATLAB: %1.10f seconds\n',tb);
fprintf('Speed factor serial/parallel: %1.10f seconds\n',ta/tap);
fprintf('Speed factor matlab/parallel: %1.10f seconds\n',tb/tap);