function f=testfun(x)
%    f=0.1+0.8*x;
%    f=1.1-1.9*x; %stable eqb, but declining func
    f=sin(x*10+15)/10+x;
%    f=sin(x*150+15)/10+x;
%     f=sin(x*55+15)/35+x;
end

