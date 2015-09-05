% Test for underflow/owerflow in lnsum and logsum formalas
clear;
clc;
mex f_lnsum.c;

%% TEST f_lnsum
%  test for overflow

fprintf('Test for overflow/underflow \n');
n=3;
v=10000;

fprintf('equal x values\n');
x=ones(n,1)*v;
for j=1:n;
    fprintf('x(%d)=%f  '  ,j,x(j));
end;
fprintf('\n');

eta=0;
logsum=f_lnsum(x,eta);
fprintf('lnsum(eta=%f)=%f\n',eta,logsum);

eta=0.000001;
logsum=f_lnsum(x,eta);
fprintf('lnsum(eta=%f)=%f\n',eta,logsum);

eta=1;
logsum=f_lnsum(x,eta);
fprintf('lnsum(eta=%f)=%f\n',eta,logsum);

abs(v);
logsum=f_lnsum(x,eta);
fprintf('lnsum(eta=%f)=%f\n',eta,logsum);

fprintf('\n');
fprintf('first element is smallest \n');
x=ones(n,1)*v;
x(1)=v-2;
for j=1:n;
    fprintf('x(%d)=%f  '  ,j,x(j));
end;
fprintf('\n');

eta=0;
logsum=f_lnsum(x,eta);
fprintf('lnsum(eta=%f)=%f\n',eta,logsum);

eta=0.000001;
logsum=f_lnsum(x,eta);
fprintf('lnsum(eta=%f)=%f\n',eta,logsum);

eta=1;
logsum=f_lnsum(x,eta);
fprintf('lnsum(eta=%f)=%f\n',eta,logsum);

eta=abs(v);
logsum=f_lnsum(x,eta);
fprintf('lnsum(eta=%f)=%f\n',eta,logsum);

fprintf('\n');
fprintf('first element is smallest \n');
x=ones(n,1)*v;
x(1)=v-2;
for j=1:n;
    fprintf('x(%d)=%f  '  ,j,x(j));
end;
fprintf('\n');

eta=0;
logsum=f_lnsum(x,eta);
fprintf('lnsum(eta=%f)=%f\n',eta,logsum);

eta=0.000001;
logsum=f_lnsum(x,eta);
fprintf('lnsum(eta=%f)=%f\n',eta,logsum);

eta=1;
logsum=f_lnsum(x,eta);
fprintf('lnsum(eta=%f)=%f\n',eta,logsum);


fprintf('\n');
fprintf('first element is largest \n');
x=ones(n,1)*v;
x(1)=v+2;
for j=1:n;
    fprintf('x(%d)=%f  '  ,j,x(j));
end;
fprintf('\n');

eta=0;
logsum=f_lnsum(x,eta);
fprintf('lnsum(eta=%f)=%f\n',eta,logsum);

eta=0.000001;
logsum=f_lnsum(x,eta);
fprintf('lnsum(eta=%f)=%f\n',eta,logsum);

eta=1;
logsum=f_lnsum(x,eta);
fprintf('lnsum(eta=%f)=%f\n',eta,logsum);

eta=abs(v);
logsum=f_lnsum(x,eta);
fprintf('lnsum(eta=%f)=%f\n',eta,logsum);

