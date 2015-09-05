function res=symmetry_checker2(g1,g2)
%compare costs and values of firm1 in game structure g1 to
%        costs and values of firm2 in game structure g2

n=numel(g1);
if n~= numel(g2)
    error 'Only similar game structures can be compared!'
end

i1a=[1 2 7 8 11 16 17 ]; %columns from solution for firm 1
i1b=[1 2 5 6];
i1=[i1a size(g1(1).solution,2)+i1b]; %firm 1 in g1
i2a=[2 1 9 10 12 18 19]; %columns from solution for firm 1
i2b=[3 4 7 8];
i2=[i2a size(g2(1).solution,2)+i2b]; %firm 2 in g2

res=[];
for i=1:n; %start with end game
    ss1=[g1(i).solution g1(i).ec];
    ss1(isnan(ss1(:,1)),:)=[];
    ss1(ss1(:,13)~=1,:)=[]; %choose selected equilibria
    ss1=ss1(:,i1);

    ss2=[g2(i).solution g2(i).ec];
    ss2(isnan(ss2(:,1)),:)=[];
    ss2(ss2(:,13)~=1,:)=[]; %choose selected equilibria
    ss2=ss2(:,i2);

    ss1=sortrows(ss1,1);
    ss2=sortrows(ss2,1);
    %output
    %res=[res;max(ss1-ss2,[],1)]; %one row per layer of the game
    %res=[res;[ss1(:,[1 2]) ss1(:,[3:end])-ss2(:,[3:end])];999*ones(1,numel(i1))];
    res=[res;[ss1(:,[1 2]) 1e-10*floor(1e10*(ss1(:,[3:end])-ss2(:,[3:end])))];999*ones(1,numel(i1))];
end
stat1=sum(res(~isnan(res)));
if abs(stat1)>eps
    fprintf('The game structures are ASYMMETRIC, max diff in values is %1.3e\nRows in output are diffs in values from the bottom up to the second layer of the game\n',abs(stat1));
else
    fprintf('The game structures are SYMMETRIC up to the tolearnce %1.3e\nRows in output are diffs in values from the bottom up to the second layer of the game\n',eps);
end


%try to output column names
try
    evalin('base',['assignin(''caller'',''gcol'',gcol)']);
    evalin('base',['assignin(''caller'',''gcol_ec'',gcol_ec)']);
    lbl1={ gcol{i1a} gcol_ec{i1b} };
    lbl2={ gcol{i2a} gcol_ec{i2b} };
    for i=1:size(i1,2)
        fprintf('%-21s %-21s\n',lbl1{i},lbl2{i})
    end
catch err
end

