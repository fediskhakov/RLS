%         Fedor Iskhakov, University Technology Sidney
%         John Rust, University of Maryland
%         Bertel Schjerning, University of Copenhagen
%         Sidney, March 2012

function res = moviemaker (par,mp,sw,parameter,interval,filename,clim)
% This function makes movies by making pay-off map for a series of parameter values for given parameter within given interval

nsteps=600; %number of steps/frames in the movie
nsteps_nodisk=120;
if isempty(filename)
    SAVETODISK=0;
else
    SAVETODISK=1;
end

% do it
if sw.esr~=99
    error 'ESR strings must be activated (sw.esr must be 99) to construct the movie'
end
close all;

%% COMPILE 
%if false
if true
    fprintf('Compiling... '); 
    tic
    mex -largeArrayDims leapfrog.c -DMAXEQB=3 -DPRINTeqbloop=0 -DPRINTeqbstr=0
    tc=toc; fprintf('Compiled model in %1.10f (seconds)\n\n',tc);
end

% make sure sw is right
sw.esrstart=0; %index of the first eqstring

%make special figure
scrsz = get(0,'ScreenSize');
mfig = figure('Color',[1 1 1],'NextPlot','replacechildren','Position',[scrsz(3)*(1-1/1.5)/2 scrsz(4)*(1-1/1.5)/2 scrsz(3)/1.5 scrsz(4)/1.5]);

if SAVETODISK
    %start movie
    %vidObj = VideoWriter(filename,'Motion JPEG 2000');
    %vidObj = VideoWriter(filename,'Uncompressed AVI');
    vidObj = VideoWriter(filename,'Motion JPEG AVI');
    vidObj.FrameRate=20;
    open(vidObj);
else
    nsteps=nsteps_nodisk;
end

%% MAIN CYCLE over parameter

for parametervalue=interval(1):(interval(2)-interval(1))/(nsteps-1):interval(2)
    %recalculate changing parameter
    eval(sprintf('%s = %16.16f;',parameter,parametervalue));

    % recalculate depend parameters
    [par mp]=f_update_params(par,mp);
    
    fprintf('MovieMaker step %s=%3.3f : ',parameter,parametervalue);

    %MONOPOLY SOLUTION
    mpm=mp;
    mpm.tpm=[1 0; 1 0];
    swm.alternate=true;
    swm.analytical=false;
    swm.esr=5;
    swm.esrmax=1;
    swm.esrstart=0;
    fprintf('monopoly ');
    tic;
    [a,b,gm,c]=leapfrog(par,mpm,swm);
    ts=toc;
    fprintf('(%1.10f seconds) : ',ts);
    mon=gm(par.nC).solution(1,7)+gm(par.nC).solution(1,9);

    tic; 
    [bne, br, g, eqbstr]=leapfrog(par,mp,sw);
    ts=toc;
    fprintf('Solved model in %1.10f (seconds)\n',ts);
    eqbstr(:,isnan(eqbstr(1,:)))=[];
    eqbstr=eqbstr';
    eqbstr(:,8)=eqbstr(:,8)/(mon);

    
    delete(allchild(mfig));%clean the figure
    % plot
    if isempty(clim)
        graph.EqbstrPlot (eqbstr,[2 8],mon,mp,sw,'Pay-off map for all found equilibria',mfig);
    else
        graph.EqbstrPlot (eqbstr,[2 8],mon,mp,sw,'Pay-off map for all found equilibria',mfig,clim);
    end
    % add text annotation
    text (.95,.95,0,sprintf('%s=%3.3f',parameter,parametervalue), ...
        'Units','normalized','VerticalAlignment','Top','HorizontalAlignment','Right','FontSize',12,'FontWeight','bold');
    %draw without buffering
    drawnow
    if SAVETODISK
        %add frame to movie
        currFrame = getframe(mfig);
        writeVideo(vidObj,currFrame);    
    end
end

if SAVETODISK
    %finish movie
    close(vidObj);
end

end %function


