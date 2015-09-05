clear
clc
load gameinfo.mat
load g.mat
fprintf('MAX(ABS(C-MATLAB))\n');
fprintf('-----------------------------------------------------------------------------------------------------------------------\n');
fprintf('            v10              v11              v20               v21             p1           p2\n');
for iC=1:1:10;
JRp=gameinfo(iC).probmat;
JRv=gameinfo(iC).payoffmat;
Cg=g(iC).solution;
Cg(:,[1:2])=Cg(:,[1:2])+1;

sel=(Cg(:,1)>=iC).*(Cg(:,2)>=iC);
sel=sel+(1-sel).*(isnan(Cg(:,1))==0);
sel=(isnan(Cg(:,1))==0);
Cg=Cg(sel==1,:);
Cvp=Cg(:,[1:2 4:5 7:12]);
% Cvp=Cg(:,[1:2 4:5 9:10 7:8 11:12]);

JRvp=[JRv(:,1:8) JRp(:,6:7)];

JRvp=sortrows(JRvp,[1:2 9]);
Cvp=sortrows(Cvp,[1:2 9]);
c1c2=Cvp(:,1:2);
edges=((Cvp(:,1)==iC)+(Cvp(:,2)==iC))>=1;
checkmat=Cvp-JRvp;
checkmat1=checkmat(edges==1,:);
checkmat2=checkmat(edges==0,:);

% fprintf('MAX(ABS(C-JR))\n');
% fprintf('-----------------------------------------------------------------------------------------------------------------------\n');

fprintf('\n');
% fprintf('FULL GAME, ic=%d\n', iC);
fprintf('ic=%d    %e    %e    %e    %e   %e  %e   \n', iC, max(abs(checkmat(:,5:end))));
% fprintf('\n');

% fprintf('EDGES, ic=%d\n', iC);
% fprintf('ic1 ic2  c1              c2            v10              v11              v20               v21             p1           p2\n');
% fprintf('%d  %d  %e   %e  %e    %e    %e    %e   %e  %e   \n', max(abs(checkmat1)));
% 
% fprintf('\n');
% fprintf('INTERIOR, ic=%d\n', iC);
% fprintf('ic1 ic2  c1              c2            v10              v11              v20               v21             p1           p2\n');
% fprintf('%d  %d  %e   %e  %e    %e    %e    %e   %e  %e   \n', max(abs(checkmat2)));

end;
