% best response for firm 1

% 

        function [drp]=df(p,c1,c2,c_og, og, sigma);


	drp=zeros(2,2);

	if (og==0);

	tmp=pr(p,c_og, og, sigma);

	drp(1,1)=-(1-tmp)-(p(1)-c1)*tmp*(1-tmp)/sigma;
	drp(2,2)=-tmp-(p(2)-c2)*tmp*(1-tmp)/sigma;
	drp(1,2)=(p(1)-c1)*tmp*(1-tmp)/sigma;
	drp(2,1)=(p(2)-c2)*tmp*(1-tmp)/sigma;

        else;

	tmp1=pr1(p,c_og, og, sigma);
	tmp2=pr2(p,c_og, og, sigma);

	drp(1,1)=-(1-tmp1)-(p(1)-c1)*tmp1*(1-tmp1)/sigma;
	drp(2,2)=-(1-tmp2)-(p(2)-c2)*(1-tmp2)*tmp2/sigma;
	drp(1,2)=-(p(1)-c1)*tmp1*tmp2/sigma;
	drp(2,1)=-(p(2)-c2)*tmp1*tmp2/sigma;

        end;
