function io=iota(g,n)
io=0;
N=size(n,1);
for j=1:N
  io=io+g(j,:)*prod(n(1:j-1,:),1);
end   
end

