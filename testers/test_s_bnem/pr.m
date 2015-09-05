% pr.m: probability of choosing firm 1 as a function of price of firm 1 and 2

    function [rpr]=pr(p, c_og, og, sigma);


       if (p(1) > p(2));

         if (sigma <= 0);

             rpr=0;

         else;

	  rpr=exp(-(p(1)-p(2))/sigma); 
	  rpr=rpr/(1+rpr);

         end;

       else;

         if (sigma <= 0);

             rpr=1;

         else;

          rpr=exp(-(p(2)-p(1))/sigma);
	  rpr=1/(1+rpr);

         end;

       end;

