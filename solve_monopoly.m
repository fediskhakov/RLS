function mon=solve_monopoly(mp,par,sw)

%MONOPOLY SOLUTION
mpm=mp;
mpm.tpm=[1 0; 1 0];
swm=sw;
swm.alternate=true;
swm.analytical=true;
swm.esr=5;
swm.esrmax=1;
swm.esrstart=0;
[a,b,gm,c]=leapfrog(par,mpm,swm);
mon=gm(par.nC).solution(1,7)+gm(par.nC).solution(1,9);

end%function
