   function [sp]=backsucc(s, numeq)

   %modify for different number of equilibria at higher nodes.

%    Input numeq;  % vector of dimension (n,1) providing the maximum number
% of equilibria in each state
n=size(s,1);
sp=s;
for i=1:n;
    
    if (numeq(i) == 1);
        sp(i)=1;
    else;
        if (sp(i)<=numeq(i)-1);
            sp(i)=sp(i)+1;
            return;
        else;
            sp(i)=mod(sp(i)+1,numeq(i));
        end;
        
    end;
    
end;
   
