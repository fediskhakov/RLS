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
par.ctol=1e-12;

% Adjust mp structure - see setup for description
mp.eta=0.15;
mex -largeArrayDims leapfrog.c 

tic;
[bne, br ,g]=leapfrog(par,mp);
t=toc;
fprintf('Model solved! Runtime=%5.15f sec.\n',t);


pause

%% SECTION 2.0: SECOND ORDER BEST RESPONSE FUNCTION COLARED BY ROOTS
c1=3.5; c2=3; iC=5; 
iF=1;
figure('Name','Second order best response function - FIRM 1');
graph.br2byroots(br, g, mp.eta, c1,c2,iC,iF, 0.001);
return


brc=br(iC).br;
i= (brc(:,3) == c1 & brc(:,4) == c2);
brc=brc(i,:);

r=3;lo=.42; hi=.91;
fx=@(x) fun.invbr2byroot(brc, mp.eta, x, iF, r);
eqb=fun.findeqb(par, lo, hi, fx);
for i=1:size(eqb,1)
    line([eqb(i,4); eqb(i,4)], [0; 1],'Color', 'r', 'LineStyle', ':');
    line([eqb(i,2); eqb(i,2)], [0; 1], 'Color', 'b','LineStyle', ':');
end

return
