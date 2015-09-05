%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                   %%%%
%%%%    This static Matlab class contains MATLAB functions generate    %%%%
%%%%    simulation for the leapfrogging paper                          %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                   %%%%
%%%% BY:     Fedor Iskhakov, University Technology Sidney              %%%%
%%%%         John Rust, University of Maryland                         %%%%
%%%%         Bertel Schjerning, University of Copenhagen               %%%%
%%%%                                                                   %%%%
%%%% THIS VERSION: March 2011                                          %%%%
%%%%                                                                   %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef simul
    properties
        c0;
        c10;
        c20;
    end
    methods (Static)
        % simul.setup: set default values of simulation parameters 
        function [sp]=setup(par)
            sp.m0=1;
            sp.c0=par.cmax;
            sp.c10=par.cmax;
            sp.c20=par.cmax;
            sp.T=500;
        end
        % sequence: simulate a single sequence of equilibrium realizations of the dynamic duopoly game
        % with leapfrogging
        function [s]=sequence(g,sp, mp, par, alternate,varargin)
            %INPUTS:
            % ...
            % [optional] seed=='randseed' regenerate seed randomly on every call

            %process optinal input
            randseed=0;
            if nargin>5
                if ischar(varargin{1})
                    if strcmp(varargin{1},'randseed')
                        randseed=1;
                    else
                        error 'Wrong optional parameter, expecting `randseed`';
                    end
                end
            end
            getrand(0);% reset the seed of the getrand

            a=g(1).solution;
            cgrid=unique(a((a(:,3)==0),4)); %grid over costs
            % Initialize state variables
            iC=find(cgrid<=sp.c0, 1, 'last');    % state of the art    
            ic1=find(cgrid<=sp.c10, 1, 'last');  % firm 1   
            ic2=find(cgrid<=sp.c20, 1, 'last');  % firm 2
            if (ic1<iC) 
                error('ERROR: Initial value of sp.c10 is smaller than state of the art cost, s.c0');
            end
            if (ic2<iC) 
                error('ERROR: Initial value of sp.c20 is smaller than state of the art cost, s.c0');
            end
            
            % initialization of output structure
            s.c    =  NaN(sp.T,1);
            s.c1   =  NaN(sp.T,1);
            s.c2   =  NaN(sp.T,1);
            s.i1   =  NaN(sp.T,1);
            s.i2   =  NaN(sp.T,1);
            s.pf1  =  NaN(sp.T,1);
            s.pf2  =  NaN(sp.T,1);
            s.m  =    NaN(sp.T,1);
            s.ip1  =  NaN(sp.T,1);
            s.ip2  =  NaN(sp.T,1);
            s.t    =  NaN(sp.T,1);
            s.Tend =  0;
            m=sp.m0;
            ti=1; % Indicator for technological improvement
            for t=1:sp.T;
                if randseed
                    u_m=rand(1,1);
                else
                    u_m=getrand(1);
                end
                s.m(t,1)=m;
                s.t(t,1)=t;
                s.c(t,1)=g(iC).c;
                s.c1(t,1)=cgrid(ic1);
                s.c2(t,1)=cgrid(ic2);
                if ti>0;
                    a=g(iC).solution;
                    seleqb=(a(:,13)==1);
                    C1=a(:,4);
                    C2=a(:,5);
                    a=a(((C1==s.c1(t,1)).*(C2==s.c2(t,1)).*(seleqb))==1,:);
                    a=a(1,:); % in case a has more than 1 element, take the first one
                            % This happens when in RLS loop there are left-overs in the g structure
                            % from previous iterations of the loop (these are not overwritten for 
                            % efficiency)  This is NOT A PROBLEM                         
                    p1=a(:,11);
                    p2=a(:,12);
                end
                if alternate==1;
                    p1=p1*(m==1);
                    p2=p2*(m==2);
                end
                s.ip1(t,1)=p1;
                s.ip2(t,1)=p2;
                if randseed
                    u_ti=rand(1,1);
                    u_inv=rand(2,1);
                else
                    u_ti=getrand(2);
                    u_inv=[getrand(3);getrand(4)];
                end
                if alternate==1;
                    %alternating move
                    if (m == 1);
                      s.i1(t,1)=(u_inv(1) <= p1);
                      s.i2(t,1)=0;
                    else
                      s.i1(t,1)=0;
                      s.i2(t,1)=(u_inv(2) <= p2);
                    end;
                else
                    %simultanous move
                      s.i1(t,1)=(u_inv(1) <= p1);
                      s.i2(t,1)=(u_inv(2) <= p2);
                end
                s.pf1(t,1)=a(:,14);
                s.pf2(t,1)=a(:,15);
                
                ti=0;
                if mp.onestep==1;
                        ti=(u_ti <= fun.ipr(s.c(t,1), mp.c_tr));
                else                
                    for jC=1:(iC-1)
                        %                     jp=(iC-jC)*nC+iC;
                        ti=ti+(u_ti <= fun.ipr(s.c(t,1), mp.c_tr)*par.pti(iC,iC-jC));
                    end
                end
                
                %  BUG HERE: WRONG CODE COMMENTED OUT. 
                %  PREVIOUSLY, FIRMS COULD IMPLEMENT THE NEXT PERIOD STATE OF THE ART
                %  RATHER THIS PERIOD STATE OF THE ART SHOULD BE
                %  IMPLEMENTED (IN THE NEXT PERIOD - TIME TO BUILD)
%                 iC=iC-ti;
%                 ic1=iC*s.i1(t,1)+(1-s.i1(t,1))*ic1;
%                 ic2=iC*s.i2(t,1)+(1-s.i2(t,1))*ic2;
               
                ic1=iC*s.i1(t,1)+(1-s.i1(t,1))*ic1;
                ic2=iC*s.i2(t,1)+(1-s.i2(t,1))*ic2;
                iC=iC-ti;

                m=1+(u_m(1) <= mp.tpm(m,2));

                %  Store result in a matrix with ngp rows, one row for each stage of the
                %  game. First columt should hold date at which technology improvemen
                %  occur, second and third colmen hold date at which firm 1 and 2 invests
                
                if ((ic1==1) || (ic2==1))
                    s.Tend=s.Tend+1;
                end
                if s.Tend>3;
                    s.Tend=t;
                    break;
                elseif t==sp.T;
                    s.Tend=t;
                end
                
            end
            s.t=s.t*mp.dt;
        end
    % next functiuon here    
    end
end

function r=getrand(column)
    persistent rand_current_index randstream1;
    if column==0
        %initialize
        rand_current_index=1;
        r=NaN;

        %initialize
        if exist('randstream.mat','file')
            %use first variable in randstream.mat
            d1=load('randstream.mat');
            d2=fieldnames(d1);
            randstream1=getfield(d1,d2{1});
            fprintf('Using stream of random number from randstream1.mat.\n');
        else
            randstream1=rand(1000,1);
            save randstream.mat randstream1;
            fprintf('Generated new randstream1 and saved it in randstream1.mat.\n');
        end
        return;
    end

    %return value
    s=numel(randstream1);
    r=randstream1(mod( (rand_current_index-1)*5+column,s));
    rand_current_index=rand_current_index+1;
end


