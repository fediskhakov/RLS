function res=lexistring(index,modul)
%if index is column-vector:
%   returns the lexicographical string when letters are 0 to modul-1
%   equivalently, returns digits of index in modul-modular arithmetics
%if index is row lexistring:
%   returns the index
%   (treats the last element as lowest digit!!!)
if isempty(index) || isempty(modul)
    res=NaN;
    return;
end
if size(index,2)==1
    %return string for index
    for i=1:numel(index)
        if index(i)==0
            res0{i}=0;
        else
            %%%%%%%% MAIN BODY
            m=[1];
            while floor(index(i)/m(1))>0
                m=[m*modul 1];
            end
            res0{i}=mod(index(i),m(1:end-1));
            res0{i}=floor(res0{i}./m(2:end));
            %%%%%%%%%%%%%%%%%%
        end
    end
    %collect results in one table
    if numel(index)==1
        res=res0{1};
    else
        res=zeros(numel(index),max(cellfun(@(x) size(x,2),res0)));
        for i=1:numel(index)
            res(i,end-numel(res0{i})+1:end)=res0{i};
        end
    end

else
    %return index for strings
    for i=1:size(index,1)
        str=index(i,:);
        res(i)=0;
        if sum(str>modul-1)>0
            warning ('Found digits which are too high for given module, skipping!');
            res(i)=NaN;
            continue;
        end
        for j=1:numel(str)
            res(i)=res(i)+str(numel(str)-j+1)*modul^(j-1);
        end
    end
    res=res';
end    
end%of function
