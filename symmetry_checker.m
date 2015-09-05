% symmetry_checker.m: checks the payoff matrices in a gameinfo structure to see if the
%                     selected equilibrium payoffs are symmetric
%
%                     John Rust, Georgetown University, April, 2012
%                     Fedor Iskhakov

function res=symmetry_checker(g)
n=numel(g);
res.v=[];
res.ec=[];
for i=1:n-1; %start with end game
    ss=g(i).solution;
    ec=g(i).ec;
    ec(isnan(ss(:,1)),:)=[];
    ss(isnan(ss(:,1)),:)=[];
    ec(ss(:,13)~=1,:)=[];%choose selected equilibria
    ss(ss(:,13)~=1,:)=[];%choose selected equilibria
    %cost my, other, my 4 value functions
    ss1=ss(:,[1 2 7 8 11 16 17]); %firm1
    ss2=ss(:,[2 1 9 10 12 18 19]); %firm2
    ec1=ec(:,[1 2 5 6]); %firm1
    ec2=ec(:,[3 4 7 8]); %firm2
    [ss1 index]=sortrows(ss1,1);
    ec1=[ss1(:,1:2) ec1(index,:)];
    [ss2 index]=sortrows(ss2,1);
    ec2=[ss2(:,1:2) ec2(index,:)];
    %output
    res.v=[res.v;max(ss1-ss2)];
    res.ec=[res.ec;max(ec1-ec2)];
end
stat1=sum(res.v(~isnan(res.v)));
stat2=sum(res.ec(~isnan(res.ec)));
if abs(stat1)>eps
    fprintf('The solution is ASYMMETRIC, max diff in values is %1.3e\nRows in output are diffs in values from the bottom up to the second layer of the game\n',abs(stat1));
else
    fprintf('The solution is SYMMETRIC up to the tolearnce %1.3e\nRows in output are diffs in values from the bottom up to the second layer of the game\n',eps);
end

if abs(stat2)>eps
    fprintf('The costs are ASYMMETRIC, max diff in values is %1.3e\nRows in output are diffs in values from the bottom up to the second layer of the game\n',abs(stat2));
else
    fprintf('The costs are SYMMETRIC up to the tolearnce %1.3e\nRows in output are diffs in values from the bottom up to the second layer of the game\n',eps);
end


%{
   payoffmat=ginfo(i).payoffmat;

   payoffmat=payoffmat(find(payoffmat(:,10)==1),3:8);

   sp=size(payoffmat,1); 

   checked=[];

   for j=1:sp;
   
     ct=find(payoffmat(:,1) == payoffmat(j,1) & payoffmat(:,2) == payoffmat(j,2));
     ct1=find(payoffmat(:,1) == payoffmat(j,2) & payoffmat(:,2) == payoffmat(j,1));

     if (payoffmat(j,1)== payoffmat(j,2));


         ctol=max([abs(payoffmat(j,3)-payoffmat(j,5)) abs(payoffmat(j,4)-payoffmat(j,6))]);

         if (ctol > tol);
            tol=ctol;
         end;

         fprintf('diagonal element (c1,c2,c)=(%g,%g,%g) tol=%g\n',payoffmat(j,1),payoffmat(j,2),ginfo(i).c,tol);

     else;

         if (size(checked,1) > 0);

         s1=sum(checked(:,1) == payoffmat(j,1) & checked(:,2) == payoffmat(j,2));
         s2=sum(checked(:,1) == payoffmat(j,2) & checked(:,2) == payoffmat(j,1));

         else;

          s1=0;
          s2=0;

         end;

         if (s1 == 0 && s2 == 0);

         checked=[checked; [payoffmat(j,1) payoffmat(j,2)]];

         ctol=max([abs(payoffmat(ct,3)-payoffmat(ct1,5)) abs(payoffmat(ct,4)-payoffmat(ct1,6))]);

         if (ctol > tol);
            tol=ctol;
         end;

         ctol=max([abs(payoffmat(ct1,3)-payoffmat(ct,5)) abs(payoffmat(ct1,4)-payoffmat(ct,6))]);

         if (ctol > tol);
            tol=ctol;
         end;
         
         fprintf('off-diagonal element (c1,c2,c)=(%g,%g,%g) rows=%i %i tol=%g\n',payoffmat(j,1),payoffmat(j,2),ginfo(i).c,ct,ct1,tol);

         end;
  
     end;

   end;

end;

%}

end %function 
