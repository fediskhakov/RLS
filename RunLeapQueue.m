% This function is used to populate the queue of parameters to run the model with in automatic fashion,
% and another part is used for running the queue.
% This script can be started many times on the same computer to run on several processors, it handles
% the queue management through the file system and can run independently

function RunLeapQueue(varargin)

if nargin==0 || ~ismember(varargin{1},{'run','populate'})
    error 'Wrong argument! Need ''run'' or ''populate'''
else
    do=varargin{1};
end
%do='populate';
%do='run';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  RUNNING THE QUEUE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(do,'run')

while true
    job=false;
    load ('queue.mat');
    for i=1:numel(queue)
        if ~queue{i}.started
            queue{i}.started=true;
            try
                save queue_tmp.mat queue;
                delete('queue.mat');
                movefile ('queue_tmp.mat','queue.mat');
            catch
                %could not delete maybe because other process is doing the same
                error "Could not get another job due to disk system error.."
            end
            mp=queue{i}.mp;
            par=queue{i}.par;
            sw=queue{i}.sw;
            fl=queue{i}.cflags;
            labl=queue{i}.label;
            clear queue;
            job=true;
            break
        end
    end
    if job
        %run the job
        RunLeap(mp,par,sw,fl,labl);
    else
        break
    end
end
fprintf('Completed the queue!\n');

end %run
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  POPULATING THE QUEUE -- has to be run manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(do,'populate')
if ~exist('queue','var')
    queue={};
end

setup;

for nC=[5]
    for eta=[0]
        for tr=[-.25 .1 .25 .5]
            for step=0:1
    
                par.nC=nC;  
                mp.eta=eta;
                mp.c_tr=tr;
                mp.onestep=step;
                
                labl=sprintf('nC=%d eta=%1.4f tr=%1.4f onestep=%d',nC,eta,tr,step);
                
                [par mp]=f_update_params(par,mp);
                
                sw.alternate=false; %alternate move or simultanious move game
                sw.analytical=false;
                sw.esr=99; %eqbstrings
                sw.esrmax=100000000; %max number of distinct equilibria to be outputed (set equal to 192736405 for n=5)
                sw.esrstart=0; %index of the first eqstring

                queue{end+1}.started=false;
                queue{end}.mp=mp;
                queue{end}.par=par;
                queue{end}.sw=sw;
                queue{end}.cflags='-DMAXEQB=3 -DPRINTeqbloop=0 -DPRINTeqbstr=0';
                %label for the run
                queue{end}.label=labl;

            end
        end
    end
end

save queue.mat queue;
clear queue;
end %populate


end%function
