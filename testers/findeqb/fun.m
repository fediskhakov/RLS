%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                   %%%%
%%%%    This static Matlab class contains the MATLAB functions         %%%%
%%%%    fun.ipr(c, c_tr): probability of technological improvement     %%%%
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

classdef fun
    methods (Static)
        % fun.ipr: probablity that production costs under the state of the art technology
        % decrease in the current period
        function ip=ipr(c, c_tr)
            ip=c_tr*c;
            ip=ip/(1.0+ip);
        end
        function [D,a0,a1]=D(a0,a1,a2, B1, eta, p);
            a0=a0-eta*B1*log(p)+eta*log(p/(1-p));  %% here pvec is P1
            a1=a1+eta*B1*log(p);
            D=a1^2-4*a0*a2;
        end
        function [D,r1,r2]=roots(a0,a1,a2, B1, eta, p);
            a0=a0-eta*B1*log(p)+eta*log(p/(1-p));  %% here pvec is P1
            a1=a1+eta*B1*log(p);
            D=a1^2-4*a0*a2;
            r1 =-(1/(2*a2))*(a1-sqrt(D));
            r2 =-(1/(2*a2))*(a1+sqrt(D));
        end
        function [pimin, pimax]=FindMinMaxBr(brc, eta, iF);
            % This function compute the minimum and the maximum of the domain of the best
            % response function for firm iF using a bracketing algorithm.
            % The program searches over the best responses, pi using a braketing algorithm that
            % updates the upper and lower bound on both the minimum and the maxmium value of pi until
            % minhi-minlo<ctol and maxhi-maxlo<ctol.
            %   maxlo: lower bound of maximum of pi
            %   maxhi: upper bound of maximum of pi
            %   minlo: lower bound of minimum of pi
            %   maxhi: upper bound of maximum of pi
            
            % Unzip br
            A0=brc(5+iF-1);
            A1=brc(7+iF-1);
            B1=brc(9);
            C0=brc(10+iF-1);
            C1=brc(12+iF-1);
            a0=brc(14+iF-1);
            a1=brc(16+iF-1);
            a2=brc(18+iF-1);
            
            maxiter=100;                       % should be included in par sturture
            ctol=0.001;                  % should be included in par sturture
            
            % Iinitialize bounds on best response
            maxhi=0.999999999;
            minlo=0.000000001;
            
            % Search for max at maxhi
            % If there are real roots at pj=maxhi the maximum is larger
            % than maxhi and there is no need for further search.
            Dmax=fun.D(a0,a1,a2, B1, eta, maxhi);
            if Dmax>0;
                maxlo=maxhi;
                fprintf('Max found at Pi=1\n');
            end
            
            
            % Find the best response of firm i when opponents when opponent's
            % investment porbability is unity, ie. pat1=pi(pj=1).
            % If there is a single root on the unit interval at pi(pj=1), then a minimum is found
            % and minhi=minlo=pat1
            
            pat1=1/(1+exp((a0+a1+a2)/eta));                       % pi(pj=1)
            [D,r1, r2]=fun.roots(a0,a1,a2, B1, eta, pat1);       % Roots at pat1
            nroots=((r1>=-0.000000001)*(r1<=1.000000001)) ...
                + ((r2>=-0.000000001)*(r2<=1.000000001));   % Nuumber of roots at the unit interval
            
            maxlo=pat1;  % The maximum can never be lower pat1
            
            % If there is a single root on the unit interval when pj=1 it must be the
            % minimum value of pi and seach should stop.
            if nroots==1;
                minhi=pat1;
                minlo=minhi;
                fprintf('Min found at Pj=1\n');
                % If there are multiple roots the minimum should be lower than pat1
            else;
                minhi=pat1;
                % Checck if initial value of minlo is infact a minium
                % (i.e.at least one real root is on the unit inteerval)
                [D,r1, r2]=fun.roots(a0,a1,a2, B1, eta, minlo);
                nroots=((r1>=-0.000000001)*(r1<=1.000000001)) ...
                    + ((r2>=-0.000000001)*(r2<=1.000000001));   % Nuumber of roots at the unit interval
                
                if nroots>=1;
                    fprintf('Min found at Pi=0\n');
                    minhi=minlo;
                end;
            end
            if (maxhi-maxlo)>ctol
                for it=1:maxiter;
                    [Dh,a0h,a1h]=fun.D(a0,a1,a2, B1, eta, maxhi-(maxhi-maxlo)/2);
                    if (Dh>0);
                        maxlo=maxhi-(maxhi-maxlo)/2;
                        diff=maxhi-maxlo;
                    else
                        maxhi=maxhi-(maxhi-maxlo)/2;
                        diff=maxhi-maxlo;
                    end
                    if abs(diff)<ctol;
                        fprintf('Max found: ');
                        fprintf('it=%3.0f maxlo=%1.4f maxhi=%1.5f minlo=%1.5f minhi=%1.5f tol=%1.16f D=%1.6f\n', it, maxlo, maxhi, minlo, minhi ,diff, Dh);
                        break;
                    end
                end
            end
            
            minhi=min(minhi,maxlo);
            
            if (minhi-minlo)>ctol;
                for it=1:maxiter;
                    [D,r1, r2]=fun.roots(a0,a1,a2, B1, eta, minlo);
                    nroots=((r1>=-0.000000001)*(r1<=1.000000001)) ...
                        + ((r2>=-0.000000001)*(r2<=1.000000001));   % Nuumber of roots at the unit interval
                    
                    if (nroots==0);
                        minlo=minhi-(minhi-minlo)/2;
                        diff=minhi-minlo;
                    else
                        maxhi=minhi-(minhi-minlo)/2;
                        diff=minhi-minlo;
                    end
                    if abs(diff)<ctol;
                        fprintf('Min found: ');
                        fprintf('it=%3.0f maxlo=%1.4f maxhi=%1.5f minlo=%1.5f minhi=%1.5f tol=%1.16f D=%1.6f\n', it, maxlo, maxhi, minlo, minhi ,diff, Dh);
                        break;
                    end
                end
            end
            pimin=minhi;
            pimax=maxlo;
        end;
        function [pimin, pimax,invbr2_min, invbr2_max]=FindMinMaxBr2(brc, mp, par, pimin, pimax, iF, r);
            % This function compute the minimum and the maximum of the domain of the second ortder best
            % response function for firm iF using a bracketing algorithm.
            % The program searches over the second order best responses, pi
            % using a bisect algorithm that.
            
            % INPUTS:
            %   brc: vector of polynomial coefficients at c1, c2, c
            %   mp:  model parameters structure (this function uses mp.eta)
            %   par: parameters structure (this function uses mp.maxit and mp.ctol)
            %   pimin and pimax: minimum and maximum values of seach region (normally pimin=0+eps, pimax=1+eps)
            %   iF firm indicator (1 or 2)
            %   r: root combination number (see more in fun.invbr2byroot)
            
            % OUTPUTS:
            %   pimin, pimax: maximum and minimum of domain of second order best response function
            %   invbr2_min, invbr2_max: inverse best responses at minimum and maximum

            % SEARCH FOR MAX
            u=pimax;
            l=pimin;
            m=(l+u)/2;
            
            [invbr2_l, lsign]=fun.invbr2byroot(brc, mp.eta, l, iF, r);
            [invbr2_m, msign]=fun.invbr2byroot(brc, mp.eta, m, iF, r);
            
            if lsign==1;
                l=u;
                invbr2_max=[];
                fprintf('Max found at Pi=1\n');
            else
                for it=1:par.maxit;
                  if (msign~=1);  %% m<max
                        l=m;
                        m=(l+u)/2;
                    else
                        u=m;
                        m=(l+u)/2;
                    end
                    [invbr2_m, msign]=fun.invbr2byroot(brc, mp.eta, m, iF, r);
                    tol=abs(l-m);
                    if abs(tol)<par.ctol;
                        fprintf('Max found in interval [%g,%g] after %d iterations\n', l, u, it);
                        break;
                    end
                end
                pimax=l;

                % SEARCH FOR MIN
                u=pimax;
                l=pimin;
                m=(l+u)/2;
                
                [invbr2_max, usign]=fun.invbr2byroot(brc, mp.eta, u, iF, r);
                [invbr2_m, msign]=fun.invbr2byroot(brc, mp.eta, m, iF, r);
                
                if usign==-1;
                    u=l;
                    invbr2_min=[];
                    fprintf('Max found at Pi=1\n');
                else
                    for it=1:par.maxit;
                        if (msign==-1);   %% m<min
                            l=m;
                            m=(l+u)/2;
                        else
                            u=m;
                            m=(l+u)/2;
                        end
                        [invbr_m, msign]=fun.invbr2byroot(brc, mp.eta, m, iF, r);
                        tol=abs(l-m);
                        if abs(tol)<par.ctol;
                            fprintf('Min found in interval [%g,%g] after %d iterations\n', l, u, it);
                            break;
                        end
                    end
                    
                end
                pimin=u;
                [invbr2_min, msign]=fun.invbr2byroot(brc, mp.eta, pimin, iF, r);
            end;
        end
        function [r1,r2]=invbr(brc, eta, p, iF);
            % This function compute the inverse best reponse function for
            % firm iF
            
            B1=brc(9);
            a0=brc(14+iF-1)-eta*B1*log(p)+eta*log(p/(1-p));  %% here pvec is P1;
            a1=brc(16+iF-1)+eta*B1*log(p);
            a2=brc(18+iF-1);
            D=a1^2-4*a0*a2;
            if D>=0;
                r1 =-(1/(2*a2))*(a1-sqrt(D));
                r2 =-(1/(2*a2))*(a1+sqrt(D));
            else
                r1=nan(1,1);
                r2=r1;
            end;
        end;

        % fun.invbr2: This function computes the inverse of the second order best reponse function for firm iF
        function [p_i_invbr2]=invbr2(brc, eta, p_i, iF)
            jF=1;
            if iF==1;
                jF=2;
            end
            [p_j(1),p_j(2)]=fun.invbr(brc, eta, p_i, iF);
            [p_i_invbr2(1),p_i_invbr2(2)]=fun.invbr(brc, eta, p_j(1), jF);
            [p_i_invbr2(3),p_i_invbr2(4)]=fun.invbr(brc, eta, p_j(2), jF);
        end;
        function [invbr2,domain_id]=invbr2byroot(brc, eta, p_i, iF, r)
            % This function computes the inverse of the second order best reponse function for
            % firm iF.
            %             ir1, ir2: indicates low or high root (-1 is low, 1 is high)
            
            domain_id=1; % no real roots: p_i > p_i_max
            
            r1=[1;-1;-1;1];
            r2=[1;1;-1;-1];
            ir1=r1(r);
            ir2=r2(r);
            
            jF=1;
            if iF==1;
                jF=2;
            end
            
            invbr2=nan(1,1);
            B1=brc(9);
            a0=brc(14+iF-1)-eta*B1*log(p_i)+eta*log(p_i/(1-p_i));  
            a1=brc(16+iF-1)+eta*B1*log(p_i);
            a2=brc(18+iF-1);
            D=a1^2-4*a0*a2;
            if D>=0;
                domain_id=-1; % real roots outside unit interval,  p_i < p_i_min 
                p_j=-(1/(2*a2))*(a1+ir1*sqrt(D));
                if ((p_j)>0 & (p_j)<1)
                    domain_id=ir1; % no real roots: p_i > p_i_max
                    a0=brc(14+jF-1)-eta*B1*log(p_j)+eta*log(p_j/(1-p_j));  
                    a1=brc(16+jF-1)+eta*B1*log(p_j);
                    a2=brc(18+jF-1);
                    D=a1^2-4*a0*a2;
                    if D>=0;
                        domain_id=-ir1; % real roots outside unit interval,  p_i < p_i_min
                        invbr2=-(1/(2*a2))*(a1+ir2*sqrt(D));
                        if ((invbr2)>0 & (invbr2)<1)
                            domain_id=0; % p_i_min <= p_i =< p_i_max
                        else
                            invbr2=nan(1,1);
                        end
                        
                    end
                end
            end
        end;
        function [eqb]=FindEqb(br, par, mp, c1,c2,iC,iF);
            % This function compute equilibrium by searching for fixed point on
            % the inverse second order best reposne function. Since the
            % inverse second order best reponse function can take maximum 4
            % values for p
            
            eqb=[];
            
            brc=br(iC).br;
            i=find(brc(:,3) == c1 & brc(:,4) == c2);
            brc=brc(i,:);
            
            for r=1:4;
                [l, u,invbr2_l, invbr2_u]=fun.FindMinMaxBr2(brc, mp,par,0.000000001, 0.999999999, iF, r);
                tol=100;
                m=(l+u)/2;
                [invbr2_m, id]=fun.invbr2byroot(brc, mp.eta, m, iF, r);
                
                lsign=(l > invbr2_l);
                usign=(u > invbr2_u);
                msign=(m>  invbr2_m);
                
                if lsign~=usign;
                    for i=1:par.maxit;
                        if (msign == lsign);
                            l=m;
                            m=(l+u)/2;
                        else;
                            u=m;
                            m=(l+u)/2;
                        end;
                        %                         fprintf('%d new interval is [%g,%g]\n',i,l,u);
                        tol=abs(l-m);
                        invbr2_m=fun.invbr2byroot(brc, mp.eta, m, iF, r);
                        msign=(m > invbr2_m);
                        if tol<par.ctol;
                            eqb=[eqb; m];
                            fprintf('%d equilibrium found for root %d in interval [%g,%g]\n',i,r, l,u);
                            break
                        end
                    end
                else
                    fprintf('No equilibrium found for root %d in interval [%g,%g]\n',r, l,u);
                end
            end
        end
        function [eqb]=FindEqb2(br, par, mp, c1,c2,iC,iF);
            % This function compute equilibrium by searching for fixed point on
            % the inverse second order best reposne function. Since the
            % inverse second order best reponse function can take maximum 4
            % values for p
                        
            eqb=[];
            
            r1=[1;-1;-1;1];
            r2=[1;1;-1;-1];
            
            brc=br(iC).br;
            i=find(brc(:,3) == c1 & brc(:,4) == c2);
            brc=brc(i,:);
            
            for r=1:4;
                slope=r1(r)*r2(r)
                l=.0000000000001;
                u=0.999999999999;
                
                invbr2_l=0;
                invbr2_u=1;
                tol=100;

                m=(l+u)/2;
                [invbr2_m, signm]=fun.invbr2byroot(brc, mp.eta, m, iF, r);
                
                for i=1:par.maxit;
                        fprintf('m=%g is below domain of br2\n',m);
                    if signm==-1;
                        l=m;
                    elseif signm==1;
                        u=m;
                        fprintf('m=%g is above domain of br2\n',m);
                    else
                        fprintf('m=%g is within domain of br2\n',m);
                        if (m>  invbr2_m); % if m is above 45 degree line
                            u=m;
                            invbr2_u=invbr2_m;
                            if slope>0;
                                invbr2_l=min(invbr2_l,invbr2_u);  %% if invbr2 is increasing
                            else
                                invbr2_l=max(invbr2_l,invbr2_u);  %% if invbr2 is decreasing
                            end
                        else
                            l=m;
                            invbr2_l=invbr2_m;
                            if slope>0;
                                invbr2_u=min(invbr2_l,invbr2_u);%% if invbr2 is increasing
                            else
                                invbr2_u=max(invbr2_l,invbr2_u);%% if invbr2 is increasing
                            end
                        end
                        
                    end
                    
%                     if (l>invbr2_l)*(signm==1);
%                         %                         (l>invbr2_l); % if l is above 45 degree line
%                         fprintf('No equilibrium found for root %d in interval [%g,%g]\n',r, l,u);
%                         break
%                     end
                    
                    fprintf('%d [invbr2_m,m] is [%g,%g]\n',i,invbr2_m, m);
                    fprintf('%d [invbr2_u,u] is [%g,%g]\n',i,invbr2_u, u);
                    fprintf('%d [invbr2_l,l] is [%g,%g]\n',i,invbr2_l, l);

                    m=(l+u)/2;
                    fprintf('%d new [l,u] interval is [%g,%g]\n',i,l,u);
                    tol=abs(l-m);
                    [invbr2_m, signm]=fun.invbr2byroot(brc, mp.eta, m, iF, r);
                    
                    if tol<par.ctol;
                        eqb=[eqb; m];
                        fprintf('%d equilibrium found for root %d in interval [%g,%g]\n',i,r, l,u);
                        break
                    end
                    %                     fprintf('No equilibrium found for root %d in interval [%g,%g]\n',r, l,u);
                end
            end
        end
        function [eqb]=Bisect(brc, par, mp, iF, l,u,r);
            % This function compute equilibrium by searching for fixed point on
            % the inverse second order best reposne function. Since the
            % inverse second order best reponse function can take maximum 4
            % values for p
                        
            eqb=[];
            tol=100;
            m=(l+u)/2;
            [invbr2_m, signm]=fun.invbr2byroot(brc, mp.eta, m, iF, r);
            
            for i=1:par.maxit;
                fprintf('m=%g is below domain of br2\n',m);
                if signm==-1;
                    l=m;
                elseif signm==1;
                    u=m;
                    fprintf('m=%g is above domain of br2\n',m);
                else
                    fprintf('m=%g is within domain of br2\n',m);
                    if (m>  invbr2_m); % if m is above 45 degree line
                        u=m;
                        invbr2_u=invbr2_m;
                    else
                        l=m;
                        invbr2_l=invbr2_m;
                    end
                end
               
                m=(l+u)/2;
                fprintf('%d new [l,u] interval is [%g,%g]\n',i,l,u);
                tol=abs(l-m);
                [invbr2_m, signm]=fun.invbr2byroot(brc, mp.eta, m, iF, r);
                
                if tol<par.ctol;
                    eqb=[eqb; m];
                    fprintf('%d equilibrium found for root %d in interval [%g,%g]\n',i,r, l,u);
                    break
                end
                %                     fprintf('No equilibrium found for root %d in interval [%g,%g]\n',r, l,u);
            end
        end
        function [eqb, signm]=SearchUnstable(brc, par, mp, iF, m, r);
        % This function compute equilibrium by searching for fixed point on
            % the inverse second order best reposne function. Since the
            % inverse second order best reponse function can take maximum 4
            % values for p
                        

            
            r1=[1;-1;-1;1];
            r2=[1;1;-1;-1];
            
            tol=100;
            
            [invbr2_m, signm]=fun.invbr2byroot(brc, mp.eta, m, iF, r);
            
            fprintf('Begin succesive approximations to find unstable equiilibria\n',m);
            for i=1:par.maxit;
                [invbr2_m, signm]=fun.invbr2byroot(brc, mp.eta, m, iF, r);
                tol=abs(m-invbr2_m);
                if signm~=0;
                    fprintf('%d NO unstable equilibria found for root %d in interval [%g,%g], signm=%g\n',i,r,invbr2_m, m, signm);
                    eqb=nan(1,1);
                    break
                end
                fprintf('%d New interval [invbr2_m,m] is [%g,%g], tol=%g\n',i,invbr2_m, m, tol);
                m=invbr2_m;
                if tol<par.ctol;
                    eqb=m;
                    fprintf('%d unstable equilibrium found for root %d in interval [%g,%g], tol=%g\n',i,r,invbr2_m, m, tol);
                    break
                end
            end
        end
        function [x, domainid]=testfun(y);
            b=10*pi/2;
            a=.8/b;
            c=1;
            x=a*cos(b*y)+c*y;
            lo=0; hi=1;
            domainid=(x>hi)-(x<lo);
        end % end testfun
        function [eqbinfo]=findeqb_test(par,lo, hi, finv)
            jmax=20;
            
            % 1: Initialize eqbinfor, to include endpoints
            eqbinfo=nan(jmax,5); % cols in eqbinfo: i, eqb, l, r;
            eqbinfo(1,1:4)=[1,lo, lo, hi];
            eqbinfo(2,1:4)=[2,hi, lo, hi];
            
            J=2;
           
            % 3: start bisection algorithm
            l=lo;
            r=hi;
            
            for i=1:par.maxit;
                m=(l+r)/2;
                fprintf('search in interval [l,r] [%g, %g], m=%g \n\n',l, r,m);
                
                % search for stable equilibira
                [x1_m, retcode]=fun.SearchStable(par, m, finv, l, r);
                if retcode == 0  % if stable equilibirum is found add row to eqbinfo
                    format short g
                    J=J+1
                    eqbinfo(J,1:4)=[J,x1_m, l, r];
                    [eqbinfo(1:J,:), sortindex]=sortrows(eqbinfo(1:J,:),2);
                    j=sortindex(J);
                    
                    m=x1_m;
                end
                
                for k=1:J
                    if eqbinfo(k,2)>max(x1_m,m)
                        eqbinfo(k,3)=max([x1_m,m,eqbinfo(k,3)]); % update l if improvement
                    elseif eqbinfo(k,2)<min(x1_m,m)
                        eqbinfo(k,4)=min([x1_m,m,eqbinfo(k,4)]); % update l if improvement
                    end
                    
                    % update tolerance
                    if k==1;
                        eqbinfo(k,5)=eqbinfo(k,4)-lo;
                    elseif k==J;
                        eqbinfo(k,5)=hi-eqbinfo(k,3);
                    else
                        eqbinfo(k,5)=eqbinfo(k-1,4)-eqbinfo(k,3);
                    end

                end
                   
                eqbinfo_print=eqbinfo(1:J,:)
                [tol, j]=max(eqbinfo(1:J,5))
                if abs(tol)>par.ctol;
                    if j==1;
                        l=eqbinfo(j,3);
                        r=eqbinfo(j,4)-2*par.ctol*2;
                    elseif j==J;
                        l=eqbinfo(j,3)+2*par.ctol*2;
                        r=eqbinfo(j,4);
                    else
                        l=eqbinfo(j,3)+2*par.ctol*2;
                        r=eqbinfo(j-1,4)-2*par.ctol*2;
                    end
                else
                    break
                end
            end;
            eqbinfo=eqbinfo(1:J,:)
        end % findeqb_test
        function [eqbinfo]=findeqb(par,lo, hi, finv)
            jmax=20;
            
            % 1: Initialize eqbinfor, to include endpoints
            eqbinfo=nan(jmax,5); % cols in eqbinfo: i, eqb, l, r, stable/unstable(1/0);
            eqbinfo(1,1:4)=[1,lo, lo, hi];
            eqbinfo(2,1:4)=[2,hi, lo, hi];
            neqb=2;
            
            
     
            
            % 2: start IRS-platinium algorithm
            j=2;
            l=lo;
            r=hi;
            for i=1:par.maxit;
                if i==1;
                    m=lo+par.ctol;
                elseif i==2;
                    m=hi-par.ctol;
                else
                    m=(l+r)/2;
                end;
                
                fprintf('search in interval [l,r] [%g, %g], m=%g \n\n',l, r,m);
                
                % search for stable equilibira
                [x1_m, retcode]=fun.SearchStable(par, m, finv, l, r);
                retcode
                
 before=eqbinfo(1:neqb,:)'
                if retcode == 0  % if stable equilibirum is found add row to eqbinfo
                    neqb=neqb+1
                    eqbinfo(neqb,1:4)=[neqb,x1_m, l, r];
                    [eqbinfo(1:neqb,:), sortindex]=sortrows(eqbinfo(1:neqb,:),2);
                    j=sortindex(neqb);
                else;
                    x1_m=m;
                end;
                
                

                % update bounds
                for k=1:neqb
                    if (eqbinfo(k,2)>max(x1_m,m)) && (retcode<=0)
                        fprintf('***** l ******');
                        eqbinfo(k,3)=max([x1_m,m,eqbinfo(k,3)]); % update l if improvement
                    elseif (eqbinfo(k,2)<min(x1_m,m)) && (retcode>=0)
                        fprintf('***** r ******');
                        eqbinfo(k,4)=min([x1_m,m,eqbinfo(k,4)]); % update r if improvement
                    else
                        fprintf('***** otherwise ******');

                    end
                end
                
                                after1=eqbinfo(1:neqb,:)'
               
                if (retcode == -1) && (j==1) % below l (during left search)(j==1)
                    eqbinfo(j,3)=max(m,eqbinfo(j,3));
                elseif (retcode ==  1) &&  (j==neqb); % above r (during right search)
                    eqbinfo(j,4)=min(m,eqbinfo(j,4));
                end;
                                after2=eqbinfo(1:neqb,:)'
                % update tolerance
                for k=2:neqb
                    if k==1;
                        eqbinfo(k,5)=eqbinfo(k,4)-eqbinfo(k,3);
                    else
                        eqbinfo(k,5)=eqbinfo(k-1,4)-eqbinfo(k,3);
                    end
                end
                
                eqbinfo_print=eqbinfo(1:neqb,:);
                [tol, j]=max(eqbinfo(2:neqb,5));
                if j==1;
fprintf('***** 1 ******');
                    l=eqbinfo(j,3);
                    r=eqbinfo(j,4)-5*par.ctol*2;
                elseif j==neqb;
fprintf('***** 2 ******');
                    l=eqbinfo(j,3)+5*par.ctol*2;
                    r=eqbinfo(j,4);
                else
fprintf('***** 3 ******');
                    l=eqbinfo(j,3)+5*par.ctol*2;
                    r=eqbinfo(j-1,4)-5*par.ctol*2;
                end
                after=eqbinfo(1:neqb,:)'
                if r-l<10*par.ctol;
                    break;
                end
                
            end
            a=eqbinfo(1:neqb,:)
        end % END OF findeqb
        function [x1, retcode]=SearchStable(par, x0, finv, l, r)
            fprintf('Search for stable equilibrium\n\n');
            fprintf('  i  m        invbr2(m)     tol\n');
            for i=1:par.maxit;
                [x1, domainid] =finv(x0);
                if domainid==1;
                    x1=r+1;
                end
                if domainid==-1;
                    x1=l-1;
                end
                
                tol=abs(x0-x1);
                fprintf('%3d %1.8f %1.8f %1.8f\n',i, x0,x1, tol);
                if tol<(par.ctol/10);
                    eqb=x1;
                    retcode=0;
                    fprintf('Stable equilibrium x1=%g found after %d iterations, tol=%g\n\n',x1, i,tol);
                    break
                end
                if x1>=r;
                    fprintf('above upper bound, r=%g, x1=%g\n\n',r, x1);
                    retcode=1;
                    break
                end
                if x1<=l;
                    retcode=-1;
                    fprintf('below lower bound, l=%g x1=%g\n\n',l, x1);
                    break
                end
                x0=x1;
            end
        end %% end of search stable
        
    end %% end fun methods
end %% end fun class

