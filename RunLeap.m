function [bne, br, g, eqbstr]=RunLeap(mp,par,sw,cflags,varargin)

% This function runs the leapfrogging model independently 
% in a separate directory 'runX' and saves there all the results 
% from the run
  
% Inputs: parameters in structures
%         any other variable to be saved in result.mat
% Output: standard from leapfrog.c (also saved to disk)
%         one line is printed on the screen

% Part 1 Parameters
par.pti=f_pti(mp, par); %make sure pti matrix is updated


% Part 2 Make directory, copy files, compile
n=0;
while exist(['run' sprintf('%05d',n)],'dir')
    n=n+1;
end
dirname=['run' sprintf('%05d',n)];
mkdir(dirname);
mkdir([dirname filesep 'code']);
copyfile('*.c',['.' filesep dirname filesep ]);
copyfile('*.m',['.' filesep dirname filesep 'code' filesep ]);
%change directory
cd (dirname);

%write info file
if isunix
    fid=fopen ('_info.txt','w+');
else
 fid=fopen ('_info.txt','wt');
end
fprintf(fid,'\nMain switches (model type) for this run of leapfrog.c:\n');
fprintf(fid,'%s',evalc('disp(sw);'));
fprintf(fid,'%s\n','Parameters of the model:');
fprintf(fid,'%s',evalc('disp(mp);disp(par);'));
fprintf(fid,'Transtion prob for m:\n');
fprintf(fid,'%s',evalc('disp(mp.tpm);'));
fprintf(fid,'Transtion prob for technology:\n');
fprintf(fid,'%s',evalc('disp(par.pti);'));
fprintf(fid,'Compile flags:\n%s\n',cflags);
%compile
fprintf([dirname ' compiling:']);
tic
out=evalc(['mex -largeArrayDims leapfrog.c ' cflags ]);
tm=toc;
fprintf(fid,'Compiled in %1.10f sec.\n',tm);
fprintf(fid,'Compiler said:\n%s\n',out);
movefile('*.c',['.' filesep 'code' filesep ]);
fprintf('%1.3fsec',tm);

% Part 3 Run the model
fprintf(' solving:');
tic
out=evalc(['[bne_' dirname ', br_' dirname ', g_' dirname ', eqbstr_' dirname ']=leapfrog(par,mp,sw);']);
tm=toc;
fprintf('%1.3fsec',tm);
fprintf(fid,'Solved in %1.10f sec.\n',tm);
fprintf(fid,'Solver said:\n%s\n',out);
%delete NaN colums in eqbstrings and transpose
fprintf(' saving:');
tic
eval(['eqbstr_' dirname '(:,isnan(eqbstr_' dirname '(1,:)))=[];']);
eval(['eqbstr_' dirname '=eqbstr_' dirname ''';']);
eval(['mp_' dirname '=mp;']);
eval(['par_' dirname '=par;']);
eval(['sw_' dirname '=sw;']);
if nargin>4
    eval(['var1_' dirname '=varargin{1};']);
    eval(['save results.mat bne_' dirname ' br_' dirname ' g_' dirname ' eqbstr_' dirname ' mp_' dirname ' par_' dirname ' sw_' dirname ' var1_' dirname ';']);
else
    eval(['save results.mat bne_' dirname ' br_' dirname ' g_' dirname ' eqbstr_' dirname ' mp_' dirname ' par_' dirname ' sw_' dirname ';']);
end
eval(['bne=bne_' dirname ';']);
eval(['br=br_' dirname ';']);
eval(['g=g_' dirname ';']);
eval(['eqbstr=eqbstr_' dirname ';']);
tm=toc;
fprintf(fid,'Results saved in %1.10f sec.\n',tm);
fprintf('%1.3fsec',tm);

%Part 4 Graphics



fprintf(' done!\n');
%close info file
fclose(fid);
%change directory back
cd ('..');
end %functions
