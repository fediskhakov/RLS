function  [ptimat]=f_pti(mp,par);        
cgrid=par.cmin:(par.cmax-par.cmin)/(par.nC-1):par.cmax;    
ptimat=zeros(par.nC,par.nC);
   for j=2:par.nC;
      if (mp.onestep);
        ptimat(j,j-1)=1;
      else;
        ptimat(j,1)=betacdf(cgrid(2)/cgrid(j),mp.beta_a,mp.beta_b);
        ptimat(j,j-1)=1-betacdf(cgrid(j-1)/cgrid(j),mp.beta_a,mp.beta_b);
        for i=2:j-2;
          ptimat(j,i)=betacdf(cgrid(i+1)/cgrid(j),mp.beta_a,mp.beta_b)-betacdf(cgrid(i)/cgrid(j),mp.beta_a,mp.beta_b);
        end;
      end;
   end;
end
