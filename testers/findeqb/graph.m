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

classdef graph
    methods (Static)
        % solution: Draw solution at stage, iC of the game
        function []=solution(g,iC, par)
            nC=par.nC;
            cmin=par.cmin;
            cmax=par.cmax;
            
            a=g(iC).solution;
            
            % set min and max for axes
            zmax=[];
            zmin=[];
            for i=1:nC-1
                zmin(i)=floor(0.1*min(min(g(i).solution(:,7:10))))*10;
                zmax(i)=floor(0.1*max(max(g(i).solution(:,7:10))))*10+10;
            end
            zmax=max(zmax);
            zmin=min(zmin);
            
            C1=a((a(:,3)==0),4);
            C2=a((a(:,3)==0),5);
            
            V10=a((a(:,13)==1),7);
            V10=reshape(V10,nC-iC+1,nC-iC+1);
            V11=a((a(:,13)==1),8);  V11=reshape(V11,nC-iC+1,nC-iC+1);
            V20=a((a(:,13)==1),9);  V20=reshape(V20,nC-iC+1,nC-iC+1);
            V21=a((a(:,13)==1),10); V21=reshape(V21,nC-iC+1,nC-iC+1);
            
            p1=a((a(:,13)==1),11);  p1=reshape(p1,nC-iC+1,nC-iC+1);
            p2=a((a(:,13)==1),12);  p2=reshape(p2,nC-iC+1,nC-iC+1);
            C1=reshape(C1,nC-iC+1,nC-iC+1);
            C2=reshape(C2,nC-iC+1,nC-iC+1);
            for i=1:6
                if i==1;
                    z=V10;
                    tit='Value function, V10';
                elseif i==2;
                    z=V11;
                    tit='Value function, V11';
                elseif i==3;
                    z=p1;
                    tit='Investment probability, firm 1';
                elseif i==4;
                    z=V20;
                    tit='Value function, V20';
                elseif i==5;
                    z=V21;
                    tit='Value function, V21';
                elseif i==6;
                    z=p2;
                    tit='Investment probability, firm 2';
                end
                subplot(2,3,i), surf(C1,C2,z);
                title(tit);
                xlabel('Marginal cost, firm 1, c1');
                ylabel('Marginal cost, firm 2, c2');
                if (i==3) || (i==6);
                    axis([cmin cmax cmin cmax 0 1]);
                else
                    axis([cmin cmax cmin cmax zmin zmax]);
                end;
            end
            drawnow;
        end
        % edges: Draw edges of solution at stage, iC of the game
        function []=edges(g,iC, par)
            nC=par.nC;
            cmin=par.cmin;
            cmax=par.cmax;
            
            a=g(iC).solution;
            
            % set min and max for axes
            ymax=[];
            ymin=[];
            for i=1:nC-1
                ymin(i)=floor(0.1*min(min(g(i).solution(:,7:10))))*10;
                ymax(i)=floor(0.1*max(max(g(i).solution(:,7:10))))*10+10;
            end
            ymax=max(ymax);
            ymin=min(ymin);
            
            a=g(iC).solution;
            C1=a((a(:,3)==0),4);
            C1=reshape(C1,nC-iC+1,nC-iC+1);
            
            
            V_c1cc=a((a(:,2)==iC-1),7:10);
            Vcc2c=a((a(:,1)==iC-1),7:10);
            subplot (1,3,1), plot(C1(1,:)',V_c1cc);
            legend('v10(c1,0,0)', 'v11(c1,0,0)', 'v20(c1,0,0)', 'v21(c1,0,0)');
            axis([cmin cmax ymin ymax]);
            xlabel('Cost of firm 1, c1');
            ylabel('Value function')
            title('Endgame value function');
            
            subplot (1,3,2), plot(C1(1,:)',Vcc2c);
            legend('v10(0,c2,0)', 'v11(0,c2,0)', 'v20(0,c2,0)', 'v11(0,c2,0)');
            axis([cmin cmax ymin ymax]);
            xlabel('Cost of firm 2, c2');
            ylabel('Value function')
            title('Endgame value function');
            
            a=g(iC).solution;
            p=a((a(:,2)==iC-1),11:12);
            subplot (1,3,3), plot(C1(1,:)',p);
            legend('p1(c1,0,0)', 'p2(c1,0,0)');
            axis([cmin cmax 0 1]);
            xlabel('Cost of firm 1, c1');
            ylabel('Value function')
            title('Investment porbabilities, p1, p2');
            
            drawnow;
        end
        function []=br(br, g, eta, c1,c2, iC, spacing);
            brc=br(iC).br;
            c=g(iC).c;
            i=find(brc(:,3) == c1 & brc(:,4) == c2);
            brc=brc(i,:);
            j=[1 1];
            pvec=0.0001:spacing:1;
            np=length(pvec);
            
            p_i_vec=nan(2*np,2);
            p_i_invbr_vec=nan(2*np,2);
            
            for p_i=pvec;
                for iF=1:2;
                    [p_i_invbr(1), p_i_invbr(2)]=fun.invbr(brc, eta, p_i, iF);
                    r=find((p_i_invbr~=nan(1,1)) & (p_i_invbr>=0) & (p_i_invbr<=1));
                    for ir=1:length(r)
                        p_i_vec(j(iF), iF)=p_i;
                        p_i_invbr_vec(j(iF), iF)=p_i_invbr(r(ir));
                        j(iF)=j(iF)+1;
                    end
                end
            end
            plot(p_i_invbr_vec(1:j(1)-1,1),p_i_vec(1:j(1)-1,1),'*r','Linewidth',1,'DisplayName','P_1(P_2), Firm 1 coefs');  % inv(P1(P2)
            hold on
            plot(p_i_vec(1:j(2)-1,2),p_i_invbr_vec(1:j(2)-1,2),'*b','Linewidth',1,'DisplayName','P_2(P_1), Firm 2 coefs');  % P2(P1)
   
            title(sprintf('Best Response Functions for Firm 1 and Firm 2 at (c_1,c_2,c)=(%3.1f,%3.1f,%3.1f) and eta=%3.2f',c1,c2,c, eta));
            xlabel('Probability Firm 2 Invests');
            ylabel('Probability Firm 1 Invests');
            axis([0 1 0 1]);
            legend('Location','East');
                        
            hold off;
        end
        function []=br2(br, g, eta, c1,c2,iC,iF, spacing);
            brc=br(iC).br;
            c=g(iC).c;
            i=find(brc(:,3) == c1 & brc(:,4) == c2);
            brc=brc(i,:);
            [pmin, pmax]=fun.FindMinMaxBr(brc, eta, iF);
            j=1;
            pvec=pmin:spacing:pmax;
            np=length(pvec);
            p_i_vec=nan(4*np,1);
            p_i_invbr2_vec=nan(4*np,1);
            for p_i=pvec;
                p_i_invbr2=fun.invbr2(brc, eta, p_i, iF);
                r=find((imag(p_i_invbr2)==0) & (p_i_invbr2>=0) & (p_i_invbr2<=1));
                for ir=1:length(r)
                    p_i_vec(j)=p_i;
                    p_i_invbr2_vec(j)=p_i_invbr2(r(ir));
                    j=j+1;
                end
            end
            
            plot(p_i_invbr2_vec(1:j-1),p_i_vec(1:j-1),'*r','Linewidth',1);
            hold on
            plot(0:0.1:1,0:0.1:1,'--m','Linewidth',2);
            title(sprintf('Second Order Best Response Function for Firm %d (c_1,c_2,c)=(%3.1f,%3.1f,%3.1f) and eta=%3.2f',iF, c1,c2,c, eta));
            xlabel(sprintf('Probability Firm %d Invests', iF));
            ylabel(sprintf('Probability Firm %d Invests', iF));
            axis([0 1 0 1]);
            hold off;
        end
        function []=br2byroots(br, g, eta, c1,c2,iC,iF, spacing);
            r1=[1;-1;-1;1];
            r2=[1;1;-1;-1];
 
            brc=br(iC).br;
            c=g(iC).c;
            i=find(brc(:,3) == c1 & brc(:,4) == c2);
            brc=brc(i,:);
            [pmin, pmax]=fun.FindMinMaxBr(brc, eta, iF);
            j=1;
            pvec=pmin:spacing:pmax;
            np=length(pvec);
            invbr2_vec=nan(np,4);
             
            for r=1:4;
                for iP=1:np;
                    invbr2(iP,:)=fun.invbr2byroot(brc, eta, pvec(iP), iF, r);
                end
%                r=find((imag(invbr2(:,ir))==0) & (invbr2(:,ir)>=0) & (invbr2(:,ir)<=1));
                i=find(isnan(invbr2)==0);
                if r==1;
                    pi1=invbr2(i);
                    pj1=pvec(i);
                elseif r==2;
                    pi2=invbr2(i);
                    pj2=pvec(i);
                elseif r==3;
                    pi3=invbr2(i);
                    pj3=pvec(i);
                elseif r==4;
                    pi4=invbr2(i);
                    pj4=pvec(i);
                end
            end
                
            
            plot(pi1,pj1,'-*r','Linewidth',1);
            hold on
            plot(pi2,pj2,'-*b','Linewidth',1);
            plot(pi3,pj3,'-*k','Linewidth',1);
            plot(pi4,pj4,'-*m','Linewidth',1);
            plot(0:0.1:1,0:0.1:1,'--g','Linewidth',2);
            title(sprintf('Second Order Best Response Function for Firm %d (c_1,c_2,c)=(%3.1f,%3.1f,%3.1f) and eta=%3.2f',iF, c1,c2,c, eta));
            xlabel(sprintf('Probability Firm %d Invests', iF));
            ylabel(sprintf('Probability Firm %d Invests', iF));
            axis([0 1 0 1]);
            legend(sprintf('[%d,%d]',r1(1),r2(1)),sprintf('[%d,%d]',r1(2),r2(2)), ...
                sprintf('[%d,%d]',r1(3),r2(3)),sprintf('[%d,%d]',r1(4),r2(4)), ...
                '45 degree line','Location','SouthEast');
            hold off;
        end 
        function []=CostSequence(s, mp, tit)
            fig=figure('name','Realized sequence: costs');
            newplot;
            hold on;
            stairs(s.t,s.c1,'r','LineWidth',3);
            stairs(s.t,s.c2,'c','LineWidth',3);
            stairs(s.t,s.c,'--k','LineWidth',2);
            title({tit sprintf('sigma=%g, eta=%g, k1=%g, k2=%g, beta=%g',mp.sigma, mp.eta, mp.k1, mp.k2, mp.df)});
            xlabel('Time');
            ylabel('Marginal Costs, Prices');
            legend('c_1','c_2','c','Location','Best');
            hold off;
        end
        function []=ProfitSequence(s, mp, tit)
            fig=figure('name','Realized sequence: Prices');
            newplot;
            hold on;
            stairs(s.t,s.pf1,'r','LineWidth',3);
            stairs(s.t,s.pf2,'c','LineWidth',3);
            stairs(s.t,s.c,'--k','LineWidth',2);
            title({tit sprintf('sigma=%g, eta=%g, k1=%g, k2=%g, beta=%g',mp.sigma, mp.eta, mp.k1, mp.k2, mp.df)});
            xlabel('Time');
            ylabel('Marginal Costs, Profits');
            legend('pi_1','pi_2','c','Location','Best');
            hold off;
        end
    end
end