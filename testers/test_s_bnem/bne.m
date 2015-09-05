
function [eqinfo]=bne(c1,c2,og, c_og,sigma);

   maxit=100;

  if (sigma > 0);

   p=max([c1 c2]');
   p=[p p]';

   pp=zeros(2,1);

   tol=.00001;
   cp=2*tol;

   i=0;
%    if (c1==0)
%        fprintf('BNE: Iteration info\n');
%        fprintf('c: %10.6f %10.6f %10.6f\n', c1, c2, c_og);
%    end
   while (cp > tol && i < maxit);

     vdf=df(p,c1,c2,c_og, og, sigma);
     vf=f(p,c1,c2, c_og, og, sigma);
     s(1)=pr(p, c_og, og, sigma);
     s(2)=1-s(1);
     pp=p-inv(df(p,c1,c2,c_og, og, sigma))*f(p,c1,c2, c_og, og, sigma);
     cp=max(abs(f(pp,c1,c2, c_og, og, sigma)));
     i=i+1; 

%      if (c1==0); 
%         
% 
%          fprintf('it: %d pi: %10.6f %10.6f ', i, s(1), s(2));
%          fprintf('x: %10.6f %10.6f ', p(1), p(2));
%          fprintf('det = %10.6f ', det(vdf));
%          fprintf(' b: %10.6f %10.6f cp=%5.2e \n ', vf(1), vf(2), cp);
%      end
     p=pp;


   end;

   s=zeros(2,1);

   if (og);

    s(1)=pr1(p,c_og, og, sigma);
    s(2)=pr2(p, c_og, og, sigma);

   else;

   s(1)=pr(p, c_og, og, sigma);
   s(2)=1-s(1);

   end;

   pf=zeros(2,1);
   pf(1)=(p(1)-c1)*s(1);
   pf(2)=(p(2)-c2)*s(2);

   eqinfo=[pf' p' s']';

 else;

   pf=zeros(2,1);
   p=max([c1 c2]');
   p=[p p]';

   if (og);

     if (c1 < c2);

        if (c1 >c_og);

          pf(1)=0;
          pf(2)=0;
          s=[0 0]';

        else;

          if (c2 < c_og);
    
            pf(1)=(c2-c1);
            pf(2)=0;
            s=[1 0]';

          else;

            pf(1)=(c_og-c1);
            pf(2)=0;
            s=[1 0]';

          end;

        end;

     else;

        if (c2 > c_og);

          pf(1)=0;
          pf(2)=0;
          s=[0 0]';

        else;

          if (c1 < c_og);
    
            pf(1)=0;
            pf(2)=(c1-c2);
            s=[0 1]';

          else;

            pf(1)=0;
            pf(2)=(c_og-c2);
            s=[0 1]';

          end;

        end;

      end;

   else;

      if (c1 < c2);

        pf(1)=(c2-c1);
        pf(2)=0;
        s=[1 0]';

      else;

        pf(1)=0;
        pf(2)=(c1-c2);
        s=[0 1]';

      end;

   end;

      eqinfo=[pf' p' s']';


 end; 
 
%    eqinfo=[pf' p' s']';
%   fprintf('sig                 %g\n',sigma);
%   fprintf('Costs               c1=%g c2=%g\n',c1,c2);
%   fprintf('Equilibrium prices: p1=%g p2=%g\n',p(1),p(2));
%   fprintf('Market Shares:      P1=%g P2=%g outside good: %g\n',s(1),s(2),1.0-s(1)-s(2));
%   fprintf('Profits             pf1=%g pf2=%g\n',pf(1),pf(2));

