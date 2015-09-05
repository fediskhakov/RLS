clear;

sigma=0.0001;
c1=0;
c2=2;
c_og=5;
og=0;
eqinfo=bne(c1,c2,og, c_og,sigma);

p=eqinfo(3:4);
vd=df(p,c1,c2,c_og, og, sigma)
vf=f(p,c1,c2, c_og, og, sigma)