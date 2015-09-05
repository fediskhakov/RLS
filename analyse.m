% this function analyses the eqbstring output of the leapfog.c and makes interesting tabulations
% INPUT: eqbstr matrix
%        [optional] 'quiet' to suppress verbose output
%		 [optional] 'notab' to suppress tabulations
%		 [optional] 'alternate' to add symmetric points in eqbstr

function res=analyse(eqbstr,varargin)
	% checking optional input parameters
	doprint=true;
	print1=true;
	if nargin>1
		for i=1:nargin-1
			if ismember(varargin{i},{'alt','alternate','duplicate','dupl'})
				eqbstr=[eqbstr;eqbstr(:,[1 2 4 3 5:end])];
				fprintf('Stats include symmetric points for alternating moves game\n');
			else
				if ismember(varargin{i},{'quiet','q','silent','qt'})
					doprint=0;
					print1=0;
				elseif ismember(varargin{i},{'notab','ntab','nt'})
					doprint=0;
				end
			end
		end			
	end

	eqblabels={'lexindex','number of repetitions','value of firm 1','value of firm 2',...
    			'pure strategy equilibria',... %5
    			'symmetric equilibria',...     %6
    			'leapfrogging equilibria',...  %7
    			'efficiency',...			   %8
    			'underinvestment'};			   %9
	% 1 OVERALL COUNTS
    res.labels={'total number of equilibria',...
    			'number of distinct pay-offs',...
    			eqblabels{[5:7 9]},...
    			'min efficiency among all',...
    			'max efficiency among all',...
    			'max efficiency attained'}';

    res.statistics = [sum(eqbstr(:,2)) size(unique(eqbstr(:,3:4),'rows'),1) eqbstr(:,2)'*eqbstr(:,[5:7 9]) ... %1 2 - 3 4 5 6
    				  min(eqbstr(:,8)) max(eqbstr(:,8)) eqbstr(:,2)'*(abs(eqbstr(:,8)-max(eqbstr(:,8)))<eps)];	   %7 8 9

    if print1
	    fprintf('%35s : %5d\n','Total number of equilibria',res.statistics(1));
	    fprintf('%35s : %5d\n','Distinct equilibria pay-offs',res.statistics(2));
	    fprintf('%35s : %5d (%6.3f%%)\n','Number of pure strategy equilibria',res.statistics(3),res.statistics(3)*100/res.statistics(1));
	    fprintf('%35s : %5d (%6.3f%%)\n','Number of symmetric equilibria',res.statistics(4),res.statistics(4)*100/res.statistics(1));
	    fprintf('%35s : %5d (%6.3f%%)\n','Number of leapfrogging equilibria',res.statistics(5),res.statistics(5)*100/res.statistics(1));
	    fprintf('%35s : %5d (%6.3f%%)\n','Equilibria with underinvestment',res.statistics(6),res.statistics(6)*100/res.statistics(1));
	    fprintf('%35s : %1.8f\n','Maximum efficiency level',res.statistics(8));
	    fprintf('%35s : %5d (%6.3f%%)\n','Equilibria with maximum efficiency',res.statistics(9),res.statistics(9)*100/res.statistics(1));
    end

    %make some tabulations
    %eqbstr cols:
	%	1 lexindex
	%	2 number of repetitions
	%	3 value firm 1
	%	4 value firm 2
	%	5 pure
	%	6 symmetric
	%	7 leapfrogging
	%	8 efficiency (maybe devided by monopoly profit)
	%   9 underinvestment

	tot=sum(eqbstr(:,2));
	%add some columns and labels
	d0=eqbstr;
	d0labels=eqblabels;
	% 10 efficiency = max
	d0=[d0 abs(d0(:,8)-max(d0(:,8)))<eps];
	d0labels={d0labels{:}, 'max efficiency'};
	% 11 efficiency = min
	%d0=[d0 abs(d0(:,8)-min(d0(:,8)))<eps];
	%d0labels={d0labels{:}, 'min efficiency'};

	%TABULATIONS to make
	% indxs=[5 6;
	% 	   5 7;
	% 	   6 7];
	%indxs=nchoosek([5 6 7 9 10],2); %ALL pairs
	indxs=nchoosek([5 6 7 9 10],2); %ALL pairs

	for k=1:size(indxs,1)
		indx=indxs(k,:);
		d=groupdata(d0,indx,[2]);
		d=d(:,[2 indx]);
    	str=['Tabulation' num2str(k) ' ' d0labels{indx(1)} ' <--> ' d0labels{indx(2)}];
		res.tab(k).label=str;
		res.tab(k).data=d;

	    if doprint
	    	fprintf(['\n' str '\n']);
	    	for i=1:numel(str);fprintf('-');end;fprintf('\n');
	    	fprintf('%-5s %-5s %10s %11s %11s %11s \n',d0labels{indx(1)}(1:min(5,end)),d0labels{indx(2)}(1:min(5,end)),'counts','%total','%first','%second');
	    	for i=1:size(d,1)
	    		fprintf('%5.0f %5.0f %10.0f %10.3f%% %10.3f%% %10.3f%%\n',d(i,2),d(i,3),d(i,1),d(i,1)*100/tot,...
	    																  d(i,1)*100/sum(d(find(d(:,2)==d(i,2)),1)),d(i,1)*100/sum(d(find(d(:,3)==d(i,3)),1)));
	    	end
	    end
	end

	function res = groupdata (data, groupby, reduce)
		% This function reduces the data in rows by grouping over given colums and perfoming operation on the given columns.
		% Inputs: data - matrix with dara in rows, columns are variables, rows are observations
		%         groupby - vector of indeces of variables to be used for grouping
		%         reduce - vector of indeces of variables to be reduced using operation within each group
		% 1 find unique combinations of groupby columns in the data
		[out m n]=unique(data(:,groupby),'rows','first');
		% m is res=data(m,:); n is data=res(n,:); 
		% 2 construct result
		res=data(m,:);
		% 3 perform operation
		for i=1:size(res,1)
		    if sum((n==i))==1
		        res(i,reduce)=data((n==i),reduce);
		    else
		        res(i,reduce)=feval('sum',data((n==i),reduce));
		    end
		end
	end
end %function