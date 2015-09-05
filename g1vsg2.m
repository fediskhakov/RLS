function res=g1vsg2(g1,g2)
% This function compared two game structures

n=numel(g1);
if n~=numel(g2)
    error 'Game structures to compare should be of the same size!'
end

res=[];

for i=1:n;
    ss1=g1(i).solution;
    ss1(isnan(ss1(:,1)),:)=[];
    ss2=g2(i).solution;
    ss2(isnan(ss2(:,1)),:)=[];
    %if different sizes
    if sum(size(ss1)~=size(ss2))~=0
      res=[res;inf*ones(1,size(ss1,2))];
    else        
        %diff
        res=[res;max(ss1-ss2,[],1)];
    end
end
stat=sum(res(~isnan(res)));
if abs(stat)>eps
    fprintf('The two game structures are NOT EQUAL, max diff in values is %1.3e\nRows in output are diffs in values from the bottom up to the top layer of the game',abs(stat));
else
    fprintf('The two game structures are EQUAL up to the tolearnce %1.3e\nRows in output are diffs in values from the bottom up to the top layer of the game',eps);
end

