%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                   %%%%
%%%%    This static Matlab class contains MATLAB functions generate    %%%%
%%%%    graphs ofr the leapfrogging paper                              %%%%
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
            sp.c0=par.cmax;
            sp.c10=par.cmax;
            sp.c20=par.cmax;
            sp.T=500;
        end
        % sequence: simulate a single sequence of equilibrium realizations of the dynamic duopoly game
        % with leapfrogging
        function [s]=sequence(g,sp, mp)
            a=g(1).solution;
            cgrid=unique(a((a(:,3)==0),4));
            
            % Initialize state variables
            iC=find(cgrid<=sp.c0, 1, 'last')    % state of the art    
            ic1=find(cgrid<=sp.c10, 1, 'last')  % firm 1   
            ic2=find(cgrid<=sp.c20, 1, 'last')  % firm 2
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
            s.ip1  =  NaN(sp.T,1);
            s.ip2  =  NaN(sp.T,1);
            s.t    =  NaN(sp.T,1);
            s.Tend =  0;
         
            ti=1; % Indicator for technological improvement
            for t=1:sp.T;
                s.t(t,1)=t;
                s.c(t,1)=g(iC).c;
                s.c1(t,1)=cgrid(ic1);
                s.c2(t,1)=cgrid(ic2);
                if ti
                    a=g(iC).solution;
                    seleqb=a(:,13)==1;
                    C1=a(:,4);
                    C2=a(:,5);
                    a=a(((C1==s.c1(t,1)).*(C2==s.c2(t,1)).*(seleqb))==1,:);
                    p1=a(:,11);
                    p2=a(:,12);
                end
                s.ip1(t,1)=p1;
                s.ip2(t,1)=p2;
                u_ti=rand(1,1);
                u_inv=rand(2,1);
                s.i1(t,1)=(u_inv(1) <= p1);
                s.i2(t,1)=(u_inv(2) <= p2);
                s.pf1(t,1)=a(:,14);
                s.pf2(t,1)=a(:,15);
                
                ti=(u_ti <= fun.ipr(s.c(t,1), mp.c_tr));
                iC=iC-ti;
                ic1=iC*s.i1(t,1)+(1-s.i1(t,1))*ic1;
                ic2=iC*s.i2(t,1)+(1-s.i2(t,1))*ic2;
                
                %  Store result in a matrix with ngp rows, one row for each stage of the
                %  game. First columt should hold date at which technology improvemen
                %  occur, second and third colmen hold date at which firm 1 and 2 invests
                
                if ((ic1==1) || (ic2==1))
                    s.Tend=s.Tend+1;
                end
                if s.Tend>2;
                    s.Tend=t;
                    break;
                elseif t==sp.T;
                    s.Tend=t;
                end
                
            end
        end
    % next functiuon here    
    end
end

