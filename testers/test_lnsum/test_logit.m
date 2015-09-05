% Test for underflow/owerflow in logit and logsum formalas
mex f_logit.c

%% TEST f_logit
%  test for overflow
clear;
clc;

fprintf('Test for overflow/underflow \n');
n=3;
v=-1000;

fprintf('equal x values\n');
x=ones(n,1)*v;
eta=0;
p=f_logit(x,eta);
for j=1:n;
    fprintf('logit(x(%d)=%f,eta=%f)=%f\n',j,x(j),eta,p(j));
end;

eta=0.000001;
p=f_logit(x,eta);
for j=1:n;
    fprintf('logit(x(%d)=%f,eta=%f)=%f\n',j,x(j),eta,p(j));
end;

eta=1;
p=f_logit(x,eta);
for j=1:n;
    fprintf('logit(x(%d)=%f,eta=%f)=%f\n',j,x(j),eta,p(j));
end;

abs(v);
p=f_logit(x,eta);
for j=1:n;
    fprintf('logit(x(%d)=%f,eta=%f)=%f\n',j,x(j),eta,p(j));
end;

fprintf('\n');
fprintf('first element is smallest\n');
eta=0;
x=ones(n,1)*v;
x(1)=v-2;
p=f_logit(x,eta);
for j=1:n;
    fprintf('logit(x(%d)=%f,eta=%f)=%f\n',j,x(j),eta,p(j));
end;

eta=0.000001;
p=f_logit(x,eta);
for j=1:n;
    fprintf('logit(x(%d)=%f,eta=%f)=%f\n',j,x(j),eta,p(j));
end;

eta=1;
p=f_logit(x,eta);
for j=1:n;
    fprintf('logit(x(%d)=%f,eta=%f)=%f\n',j,x(j),eta,p(j));
end;

eta=abs(v);
p=f_logit(x,eta);
for j=1:n;
    fprintf('logit(x(%d)=%f,eta=%f)=%f\n',j,x(j),eta,p(j));
end;


fprintf('\n');
fprintf('first element is largest\n');
eta=0;
x=ones(n,1)*v;
x(1)=v+2;
p=f_logit(x,eta);
for j=1:n;
    fprintf('logit(x(%d)=%f,eta=%f)=%f\n',j,x(j),eta,p(j));
end;

eta=0.000001;
p=f_logit(x,eta);
for j=1:n;
    fprintf('logit(x(%d)=%f,eta=%f)=%f\n',j,x(j),eta,p(j));
end;

eta=1;
p=f_logit(x,eta);
for j=1:n;
    fprintf('logit(x(%d)=%f,eta=%f)=%f\n',j,x(j),eta,p(j));
end;

eta=abs(v);
p=f_logit(x,eta);
for j=1:n;
    fprintf('logit(x(%d)=%f,eta=%f)=%f\n',j,x(j),eta,p(j));
end;


