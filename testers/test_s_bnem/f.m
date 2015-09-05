% best response for firm 1

% 

        function [brp]=f(p,c1,c2, c_og, og, sigma);

	brp=zeros(2,1);

	if (og==0);

	tmp=pr(p, c_og, og, sigma);
        brp(1)=sigma-(p(1)-c1)*(1-tmp);
        brp(2)=sigma-(p(2)-c2)*tmp;

        else;

	tmp1=pr1(p, c_og, og, sigma);
	tmp2=pr2(p, c_og, og, sigma);
        brp(1)=sigma-(p(1)-c1)*(1-tmp1);
        brp(2)=sigma-(p(2)-c2)*(1-tmp2);


        end;

