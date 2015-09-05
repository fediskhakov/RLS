clc
%% COMPARE CENTRAL PARAMETRES IN C-CODE ANS MATLAB-CODE
clear 
global sw par mp dt bet k1 k2 sigma eta beta_a beta_b onestep mtp ic 
cd 'C:\Users\okobsc\Dropbox\IRS\leapfrogging\code\bertel\leapfrogging_sydney';
setup
mp.k1=1; mp.k2=10; mp.c_tr=0.05; mp.eta=0.5;  mp.eta=.0; 
mp.k1=30; mp.k2=10; mp.c_tr=1;   mp.eta=0.5;  mp.eta=.0; %% solution looks wierd and assymetric, but match John
mp.onestep=1;    
mp.tpm=[0.5 0.5; 0.5  0.5];

mp.dt=1;mp.df=exp(-0.05*mp.dt);
par.nC=11;
par.pti=f_pti(mp, par);

%%
cd 'C:\Users\okobsc\Dropbox\IRS\leapfrogging\code\bertel\john_new';
setup
checkparms(1)=(dt==mp.dt);
checkparms(2)=(bet==mp.df);
checkparms(3)=(k1==mp.k1);
checkparms(4)=(k2==mp.k2);
checkparms(5)=(sigma==mp.sigma);
checkparms(6)=(eta==mp.eta);
checkparms(7)=(beta_a==mp.beta_a);
checkparms(8)=(beta_b==mp.beta_b);
checkparms(9)=(ic==mp.c_tr);
checkparms(10)=(onestep==mp.onestep);
checkparms(11)=(sum(sum(mtp~=mp.tpm))==0);    %% NOTE MTPS ARE TRRANSPOSED IN JOHNS CODE
checkparms
%% SOLVE  MODEL USING C-CODE
clear all
cd 'C:\Users\okobsc\Dropbox\IRS\leapfrogging\code\bertel\leapfrogging_sydney';
setup;
mp.k1=5; mp.k2=10; mp.c_tr=-1; mp.eta=.0;  

mp.eta=.0; mp.onestep=1;       % MATCH    

mp.k1=1; mp.k2=10; mp.c_tr=0.05; mp.eta=0.5;  mp.eta=.0; 
mp.onestep=1;    
mp.tpm=[0.5 0.5; 0.5  0.5];


mp.k1=30; mp.k2=10; mp.c_tr=1;   mp.eta=0.5;  mp.eta=.0; %% solution looks wierd

mp.dt=1;mp.df=exp(-0.05*mp.dt);
par.nC=11;
par.pti=f_pti(mp, par);


tic
fprintf('Compiling... ')
mex -largeArrayDims leapfrog.c -DMAXEQB=3
toc
tic
nrep=1;
sw.alternate=false; %alternate move or simultanious move game
sw.esr=3; %equilibrium selection rule to be used: see setup.m or esr.c
sw.esrmax=10; %run N feasible eqstrings 
sw.esrstart=0; %index of the first eqstring
fprintf('Running... ')
for i=1:nrep;
    [bne, br, g, eqbstrings]=leapfrog(par,mp,sw);
end
ts=toc;
fprintf(' time to solve model (serial code): %1.10f (1/%g seconds)\n',ts, nrep);
save(['testers' filesep 'JRvsC_am' filesep 'g.mat'],'g');
a=g(1).solution;

save(['testers' filesep 'JRvsC_am' filesep 'g.mat'],'g');


%%
clear all
cd 'C:\Users\okobsc\Dropbox\IRS\leapfrogging\code\bertel\john_new';
regs;
%%
clear
load gameinfo0.mat  %% simultaneours move
gameinfo=gameinfo0;
%%
cd 'C:\Users\okobsc\Dropbox\IRS\leapfrogging\code\bertel\leapfrogging_sydney\testers\JRvsC_am';

load g.mat

cd 'C:\Users\okobsc\Dropbox\IRS\leapfrogging\code\bertel\leapfrogging_sydney';
symmetry_checker(g)

nC=11;
fprintf('NUBMER OF EQUILIBRIA (C)\n');
for iC=1:nC
    Cg=g(iC).solution;
    Cg(:,[1:2])=Cg(:,[1:2])+1;
    
    sel=(Cg(:,1)>=iC).*(Cg(:,2)>=iC);
    sel=sel+(1-sel).*(isnan(Cg(:,1))==0);
    sel=(isnan(Cg(:,1))==0);
    Cg=Cg(sel==1,:);
    fprintf('iC=%g,    MAX(neqb)=%g\n',iC, max(Cg(:,20)));
end

mae=[];

fprintf('MAX(ABS(C-MATLAB))\n');
fprintf('-----------------------------------------------------------------------------------------------------------------------\n');
fprintf('       v10          v11         v20          v21        x10        x11         x20        x21         p1        p2     \n');
for iC=1:nC-1;
    JRp=gameinfo(iC).probmat;
    JRv=gameinfo(iC).payoffmat;
    Cg=g(iC).solution;
    Cg(:,[1:2])=Cg(:,[1:2])+1;
    
    sel=(Cg(:,1)>=iC).*(Cg(:,2)>=iC);
    sel=sel+(1-sel).*(isnan(Cg(:,1))==0);
    sel=(isnan(Cg(:,1))==0);
    Cg=Cg(sel==1,:);
    
    Cvp=Cg(:,[1:2 4:5 7:10 11:12]);
    
    if iC==1;
        JRvp=[JRv(:,[1:4 5:8]) JRp(:,6:7)];
    else
        JRvp=[JRv(:,[1:4 5:8]) JRp(:,6:7)];
    end
    JRvp=sortrows(JRvp,[1:2 8]);
    Cvp=sortrows(Cvp,[1:2 8]);
    
    
    c1c2=Cvp(:,1:2);
    edges=((Cvp(:,1)==iC)+(Cvp(:,2)==iC))>=1;
    checkmat=[Cvp(:,1:4) Cvp-JRvp];
    checkmat1=checkmat(edges==1,[1:4 9:end]);
    checkmat2=checkmat(edges==0,[1:4 9:end]);
    JR1=JRvp(edges==1,:);
    JR2=JRvp(edges==0,:);
    C1=Cvp(edges==1,:);
    C2=Cvp(edges==0,:);
    
    % fprintf('MAX(ABS(C-JR))\n');
    % fprintf('-----------------------------------------------------------------------------------------------------------------------\n');
     fprintf('\n');   
    fprintf('\n');
    fprintf('EDGES, ic=%d\n', iC);
    fprintf('ic=%d    %1.1e    %1.1e    %1.1e    %1.1e   %1.1e  %1.1e   %1.1e    %1.1e   %1.1e   %1.1e   \n', iC, max(abs(checkmat1(:,5:end))));
%     fprintf('ic=%d    %g    %g    %g    %g   %g  %g   %g    %g   %g   %g   \n', iC, max(abs(checkmat(:,9:end))));

    
    fprintf('\n');
    fprintf('INTERIOR, ic=%d\n', iC);
    fprintf('ic=%d    %1.1e    %1.1e    %1.1e    %1.1e   %1.1e  %1.1e   %1.1e    %1.1e   %1.1e   %1.1e   \n', iC, max(abs(checkmat2(:,5:end))));
%     fprintf('ic=%d    %g    %g    %g    %g   %g  %g   %g    %g   %g   %g   \n', iC, max(abs(checkmat(:,9:end))));
    
    
    mae=[mae; max(abs(checkmat(:,9:end)))];
    
end;

fprintf('\n\nMAX(ABS(C-JR)) sumarized over all layers\n');
fprintf('-----------------------------------------------------------------------------------------------------------------------\n');
fprintf('       v10          v11         v20          v21        x10        x11         x20        x21         p1        p2     \n');
fprintf('%g    %g    %g    %g   %g  %g   %g    %g   %g   %g   \n', max(mae));
 