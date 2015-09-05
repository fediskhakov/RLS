% pr1.m: probability of choosing firm 1 as a function of price of firm 1 and 2
%        (relevant when there is an outside good, otherwise pr1=pr, where pr
%         is the probability of choosing firm 1 when there is no outside good)

    function [rpr]=pr1(p, c_og, og, sigma);

       if (p(1) > p(2));

	  rpr=exp(-(p(1)-p(2))/sigma); 

          if (p(2) < c_og);

	  rpr=rpr/(1+exp((p(2)-c_og)/sigma)+rpr);

          else;

          tmp=exp((c_og-p(2))/sigma);

	  rpr=(rpr*tmp)/(1+tmp+rpr*tmp);

          end;

       else;

          rpr=exp(-(p(2)-p(1))/sigma);

          if (p(1) < c_og);

   	    rpr=1/(1+exp((p(1)-c_og)/sigma)+rpr);

          else;

            tmp=exp((c_og-p(1))/sigma);
   	    rpr=tmp/(1+tmp+tmp*rpr);

          end;

       end;

