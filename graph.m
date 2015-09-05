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
        function []=solution(g,iC, par, plottype, fighandle)
            nC=par.nC;
            cmin=par.cmin;
            cmax=par.cmax;
            
            if plottype ==1;
                idx=[7:10 11:12];
            else
                idx=[7:8 18:19 11 12];
            end
            
            
            a=g(iC).solution;
            
            % set min and max for axes
            zmax=[];
            zmin=[];
            for i=1:nC-1
                zmin(i)=floor(0.1*min(min(g(i).solution(:,idx(1:4)))))*10;
                zmax(i)=floor(0.1*max(max(g(i).solution(:,idx(1:4)))))*10+10;
            end
            zmax=max(zmax);
            zmin=min(zmin);
            
            C1=a((a(:,3)==0),4);
            C2=a((a(:,3)==0),5);
            
            V10=a((a(:,13)==1),idx(1));
            V10=reshape(V10,nC-iC+1,nC-iC+1);
            V11=a((a(:,13)==1),idx(2));  V11=reshape(V11,nC-iC+1,nC-iC+1);
            V20=a((a(:,13)==1),idx(3));  V20=reshape(V20,nC-iC+1,nC-iC+1);
            V21=a((a(:,13)==1),idx(4)); V21=reshape(V21,nC-iC+1,nC-iC+1);
            
            p1=a((a(:,13)==1),idx(5));  p1=reshape(p1,nC-iC+1,nC-iC+1);
            p2=a((a(:,13)==1),idx(6));  p2=reshape(p2,nC-iC+1,nC-iC+1);
            
            if plottype ==2;
                p2=p2*0;
            end
            
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
                    if plottype==1
                        tit='Value function, V20';
                    else
                        tit='Value function, X20';
                    end
                elseif i==5;
                    z=V21;
                    if plottype==1
                        tit='Value function, V21';
                    else
                        tit='Value function, X21';
                    end
                elseif i==6;
                    z=p2;
                    tit='Investment probability, firm 2';
                end
                axes1=subplot(2,3,i,'Parent',fighandle);
                surf(C1,C2,z,'Parent',axes1);
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
            %             C1=a((a(:,3)==0),4);
            %             C1=reshape(C1,nC-iC+1,nC-iC+1);
            C1=a((a(:,2)==0),4);
            
            V_c1cc=a((a(:,2)==iC-1),7:10);
            Vcc2c=a((a(:,1)==iC-1),7:10);
            
            figure2 = figure('Color',[1 1 1],'NextPlot','new');
            
            axes1=subplot (1,3,1,'Parent',figure2), plot(C1',V_c1cc,'Parent',axes1);
            legend('v10(c1,0,0)', 'v11(c1,0,0)', 'v20(c1,0,0)', 'v21(c1,0,0)');
            axis([cmin cmax ymin ymax]);
            xlabel('Cost of firm 1, c1');
            ylabel('Value function')
            title('Endgame value function');
            
            axes2=subplot (1,3,2), plot(C1',Vcc2c,'Parent',axes2);
            legend('v10(0,c2,0)', 'v11(0,c2,0)', 'v20(0,c2,0)', 'v11(0,c2,0)');
            axis([cmin cmax ymin ymax]);
            xlabel('Cost of firm 2, c2');
            ylabel('Value function')
            title('Endgame value function');
            
            a=g(iC).solution;
            p=a((a(:,2)==iC-1),11:12);
            axes3=subplot (1,3,3), plot(C1',p,'Parent',axes3);
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
            if numel(brc)==0
                error('graph.br2: passed c1,c2 are not found in br().br');
            end
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
            plot(p_i_invbr_vec(1:j(1)-1,1),p_i_vec(1:j(1)-1,1),'*r','Linewidth',1);  % inv(P1(P2)
            hold on
            plot(p_i_vec(1:j(2)-1,2),p_i_invbr_vec(1:j(2)-1,2),'*b','Linewidth',1);  % P2(P1)
            
            title(sprintf('Best Response Functions for Firm 1 and Firm 2 at (c_1,c_2,c)=(%3.1f,%3.1f,%3.1f) and eta=%3.2f',c1,c2,c, eta));
            xlabel('Probability Firm 2 Invests');
            ylabel('Probability Firm 1 Invests');
            axis([0 1 0 1]);
            legend('P_1(P_2), Firm 1 coefs','P_2(P_1), Firm 2 coefs', 'Location','East');
            
            hold off;
        end

        function num_equilibria=br0eta(br,g,ic1,ic2,iC);
            % This function makes the plot of two best response
            % functions for the two firms on the unit square
            % INPUTS: solution br, g
            %         where best responses should be plotted
            %         ic1,ic2,iC - indices of the point in the state, base1 (between 1=bottom and nC=apex)
            %
            % br - solution data for the interior points

            if ic1<iC | ic2<iC
                error 'Wrong call: expecting ic1>=iC, ic2>=iC'
            end

            % compute by firms
            for iF=1:2
                PiF{iF}=[];
                PjF{iF}=[];
                if ic1==iC | ic2==iC
                    %corners and edges
                    if iC==ic1 && iC==ic2
                        %corner
                        PjF{iF}(1)=0; % x
                        PjF{iF}(2)=1; % x
                        PiF{iF}(1)=0; % y
                        PiF{iF}(2)=0; % y
                    else
                        %edge
                        if iC==1
                            %bottom layer: 0 invesment
                            PjF{iF}(1)=0; % x
                            PjF{iF}(2)=1; % x
                            PiF{iF}(1)=0;
                            PiF{iF}(2)=0;
                        else
                            %higher layers: need to compare values for one of the firms
                            if (iF==1 & ic1==iC) | (iF==2 & ic2==iC)
                                PjF{iF}(1)=0; % x
                                PjF{iF}(2)=1; % x
                                PiF{iF}(1)=0;
                                PiF{iF}(2)=0;
                            else
                                i=find(g(iC).solution(:,1) == ic1-1 & g(iC).solution(:,2) == ic2-1); %base0
                                % columns in g: 7,8 firm 1 N,I, 9,10 firm 2 N,I values
                                if (iF==1 & g(iC).solution(i,7)>g(iC).solution(i,8)) | ...
                                   (iF==2 & g(iC).solution(i,9)>g(iC).solution(i,10)) 
                                    %non investing is better
                                    PjF{iF}(1)=0; % x
                                    PjF{iF}(2)=1; % x
                                    PiF{iF}(1)=0;
                                    PiF{iF}(2)=0;
                                else
                                    PjF{iF}(1)=0; % x
                                    PjF{iF}(2)=1; % x
                                    PiF{iF}(1)=1;
                                    PiF{iF}(2)=1;
                                end
                            end
                        end
                    end
                else
                    %interior points
                    brc=br(iC).br;
                    c=g(iC).c;
                    i=find(brc(:,1) == ic1-1 & brc(:,2) == ic2-1); %base0 in br(iC).br
                    brc=brc(i,:);
                    %coefficients of the quadratic form of the best responce
                    a0=brc(14+iF-1); %constant
                    a1=brc(16+iF-1); %linear
                    a2=brc(18+iF-1); %quadratic in opponents investment probability
                    %start making the best correspondence correspondence
                    k=1;
                    PjF{iF}(k)=0; %x-axes value
                    if a2>0
                        PiF{iF}(k)=0; %with positive quadratic coef start at certain investment
                    else
                        PiF{iF}(k)=1; %otherwise with certain no investment
                    end
                    %solve for the roots for the best response correspondence
                    D=a1^2-4*a0*a2;
                    if D<0
                        %negative discriminant ==> no roots
                        %finish the line at the same value
                        PjF{iF}(k+1)=1;
                        PiF{iF}(k+1)=PiF{iF}(k);
                        k=k+1;%total points
                    else
                        %roots
                        r1 =-(1/(2*a2))*(a1-sqrt(D));
                        r2 =-(1/(2*a2))*(a1+sqrt(D));
                        %sorted roots
                        root(1)=min(r1,r2);
                        root(2)=max(r1,r2);

                        % combine the roots into the lines
                        for ir=1:2
                            if root(ir)<0
                                %root before unit interval --> replace the point from 1 to 0 and vise versa
                                PiF{iF}(k)=1-PiF{iF}(k);
                            elseif root(ir)>=0 & root(ir) <=1
                                %root inside of unit interval --> make vertical line
                                PjF{iF}(k+1:k+2)=[root(ir) root(ir)];
                                PiF{iF}(k+1:k+2)=[PiF{iF}(k) 1-PiF{iF}(k)];
                                k=k+2;
                            else
                                %root after unit interval --> ignore
                            end
                        end
                        %last point
                        PjF{iF}(k+1)=1;
                        PiF{iF}(k+1)=PiF{iF}(k);
                        k=k+1;%total points
                    end
                end
            end

            % Onto the best response function overlay the computed equilibria
            % In addition to making nice graph, this makes the visual check of correctness
            i=find(g(iC).solution(:,1)==ic1-1 & g(iC).solution(:,2)==ic2-1);
            eqbs=g(iC).solution(i,[11 12]);

            % plot
            linethick=5;
            fontsize=24;
            fig1=figure('Color',[1 1 1],'NextPlot','new');
            ax=axes('Parent',fig1);
            plot(ax,PjF{1},PiF{1},'Linewidth',linethick,'Color','red');% best response of firm 1
            hold(ax,'on');
            plot(ax,PiF{2},PjF{2},'Linewidth',linethick,'Color','black','LineStyle','--');% best response of firm 2, swapping axes
            hold(ax,'on');
            scatter(eqbs(:,2),eqbs(:,1),120,[0 0 1],'filled','Linewidth',linethick,'MarkerFaceColor','white','MarkerEdgeColor','black');
            % title(ax,'Best response functions for firm 1 (solid) and firm 2 (dashed)', ...
            %     'FontSize',fontsize);
            set(ax,'XTick',[0 .5 1],'XTickLabel',{'0','','1'}, ...
                   'YTick',[0 .5 1],'YTickLabel',{'0','','1'}, ...
                   'XLim',[-.01 1.01],'YLim',[-.01,1.01],'box','off', ...
                   'TickDir','out','DataAspectRatio',[1 1 1], ...
                   'PlotBoxAspectRatio',[1 1 1],'FontSize',fontsize,...
                   'Linewidth',linethick/2);
            xlabel(ax,'P2','FontSize',fontsize);
            ylabel(ax,'P1','FontSize',fontsize);
            hold(ax,'off');
        end

        
        function []=br2(br, g, eta, c1,c2,iC,iF, spacing);
            brc=br(iC).br;
            c=g(iC).c;
            i=find(brc(:,3) == c1 & brc(:,4) == c2);
            brc=brc(i,:);
            if numel(brc)==0
                error('graph.br2: passed c1,c2 are not found in br().br');
            end
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
        
        function []=CostSequence(s, mp, tit,varargin)
            fig=figure('name','Realized sequence: costs','NextPlot','new','Color',[1 1 1]);
            axes1=axes('Parent',fig,'FontSize',14,'FontName','Times New Roman');
 
            hold on;

            %stairs(s.t,s.c1,'r','LineWidth',3);
            %stairs(s.t,s.c2,'c','LineWidth',3);
            %stairs(s.t,s.c,'--k','LineWidth',2);

            if nargin <= 4
                pl=plot(s.t,s.c1,s.t,s.c2,'Parent',axes1,'MarkerFaceColor',[1 1 1],'LineWidth',2);
            else
                c1=s.c2;
                pl=plot(s.t,s.c1,s.t,s.c2,s.t,s.cmon,'Parent',axes1,'MarkerFaceColor',[1 1 1],'LineWidth',2);
            end

            set(pl(1),'MarkerSize',4,'Marker','o','Color',[0 0.5 0],'DisplayName','Firm 1 cost (c_1)');
            set(pl(2),'MarkerSize',3,'Marker','square','Color',[1 0 0],'DisplayName','Firm 2 cost (c_2)');
            stairs(s.t,s.c,'Parent',axes1,'LineWidth',2,'DisplayName','state of the art cost (c)','Color',[0 0 0]);
            if nargin > 4
                set(pl(3),'MarkerSize',3, 'Marker','square','Color',[0 0 1],'DisplayName','Monopolist');
                lg=legend('c_1','c_2','c_{monopoly}','c', 'Location','NorthEast');
            else
                lg=legend('c_1','c_2','c', 'Location','NorthEast');
            end

            if nargin > 3
                if strcmp(varargin(1),'');
                    title({tit});
                else
                    title({tit varargin{1}});
                end
            else
                title({tit sprintf('sigma=%g, eta=%g, k1=%g, k2=%g, dt=%g, beta=%g',mp.sigma, mp.eta, mp.k1, mp.k2,mp.dt,mp.df)});
            end
            xlabel('Time');
            ylabel('Marginal Costs, Prices');
            set(lg,'Box','off');
            ylm=get(axes1,'Ylim');
            set(axes1,'YLim',[ylm(1) ylm(2)+0.1]);
            box(axes1,'on');
            hold off;
        end
        
        function []=ProfitSequence(s, mp, tit)
            fig=figure('name','Realized sequence: Prices');
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
        
        function [statistics labels]=EqbstrCDFPlot(eqbstr,vars,dummystr,dummylbl, mp,sw,tit,xlow,varargin)
            % Plots the CDF of distribution of efficiency of equilibria in eqbstrings
            % Assumes: eqbstring(1) - lexicographical index
            %          eqbstring(2) - count of repeated equilibria
            %          eqbstring(3) - value for x axis
            %          eqbstring(4) - value for y axis
            %          eqbstring(5,6,..) - additional statistics
            % Inputs: eqbstring - as outputed by solver
            %         vars = 1xn vector of indices of the columns of eqbstring to condition on.
            %                n CDFs are plotted - one for each element in vars
            %                if vars(i)= 0 , CDF of all equilibria are shown
            %                if vars(i)< 0 , CDF of equilibria where eqbstr(:,abs(vars(i)))=1
            %                if vars(i)> 0 , CDF of equilibria where eqbstr(:,vars(i))=0
            %         sw - parameters as input for solver
            %         mp - same
            %         tit - title for the graph
            %    optional handle for the figure where graphs should be
            % Output: statistics on the found equilibria
            %         cell array of some labels for the statistics
            
            %do statistics first
            labels={'total number of equilibria';'number of distinct pay-offs';'pure strategy';'symmetric';'leapfrog';'efficiency';'underinvest'};
            title3={'','number of repetitions','value of firm 1','value of firm 2',labels{3:end}};
            labels_n={'total number of equilibria';'number of distinct pay-offs';'mixed strategy';'asymmetric';'non-leapfrogging';'efficiency'};
            title3_n={'','number of repetitions','value of firm 1','value of firm 2',labels_n{3:end}};
            
            
            n=numel(vars);
            nd=size(dummylbl,2);
            if nd>0;
                dummy=eval(dummystr);
            end         
            legendtxt='';
            for i=1:n+nd;
                if i<=n
                    if vars(i)==0;
                        mask=ones(size(eqbstr,1),1)==1;
                    else
                        if vars(i)<0
                            mask=eqbstr(:,abs(vars(i)))==0;
                        else
                            mask=eqbstr(:,vars(i))==1;
                        end
                    end
                else
                    mask=dummy(:,i-n)==1;
                end
                seqbstr=sortrows(eqbstr(mask,:),8);
                statistics(i,:) = [sum(seqbstr(:,2)) size(unique(seqbstr(:,3:4),'rows'),1) seqbstr(:,2)'*seqbstr(:,5:end)];
                
                
                if ~isempty(seqbstr)
                    Ni=[0; seqbstr(:,2)];
                    Xi=[seqbstr(1,8); seqbstr(:,8)];
                else
                    Ni=nan(1,1);
                    Xi=nan(1,1);
                end

                %plot
                if i==1;
                    if nargin>8 & ishandle(varargin{1})
                        figure1=varargin{1};
                    else
                        scrsz = get(0,'ScreenSize');
                        %figure1 = figure('Color',[1 1 1],'NextPlot','new');
                       % figure1 =figure('Color',[1 1 1],'NextPlot','replacechildren','Position',[scrsz(3)*(1-1/1.5)/2 scrsz(4)*(1-1/1.5)/2 scrsz(3)/2 scrsz(4)/2]);
                        figure1 =figure('Color',[1 1 1],'NextPlot','replacechildren');
                    end
                    axes1 = axes('Parent',figure1,'FontSize',14,'FontName','Times New Roman',...
                        'YGrid','on',...
                        'YColor',[0.25 0.25 0.25],...
                        'XGrid','on',...
                        'XColor',[0.25 0.25 0.25],...
                        'PlotBoxAspectRatio',[1 1 1]);
                    colormap('jet');
                    box(axes1,'on');
                    hold(axes1,'all');
                end
                
                % stairs(Xi, Ni)
                if i<=n;
                    if vars(i)==0
                        tmp=stairs(Xi, cumsum(Ni)/sum(Ni),'LineWidth',2);
                        set(tmp,'DisplayName', 'all equilibria');
                    else
                        if vars(i)<0
                            tmp=stairs(Xi, cumsum(Ni)/sum(Ni),'LineWidth',1);
                            set(tmp,'DisplayName', title3_n{abs(vars(i))});
                        else
                            tmp=stairs(Xi, cumsum(Ni)/sum(Ni),'LineWidth',1);
                            set(tmp,'DisplayName', title3{vars(i)});
                        end
                    end
                    
                else
                    tmp=stairs(Xi, cumsum(Ni)/sum(Ni),'LineWidth',1);
                    set(tmp,'DisplayName', dummylbl{i-n});
                end
            end
            
            xlim=get(axes1,'XLim');
            xlim(1)=min(xlow, xlim(1));
            xlim(2)=1;
            set(axes1,'XLim',xlim);
            set(axes1,'YLim',[0 1]);
            
            legend1=legend('show');
            set(legend1, 'Location','West');
            set(legend1, 'EdgeColor', [1 1 1]);
            set(legend1, 'Color', [1 1 1])
            
            title2='';
            for i=1:n+nd;
                if i<=n;
                    if vars(i)==0
                        lbl='equilibria'
                    elseif vars(i)<0
                        lbl=title3_n{abs(vars(i))};
                    else
                        lbl=title3{(vars(i))};
                    end
                else
                      lbl=sprintf('%s', dummylbl{i-n});
                end
                if i<n+nd;
                    title2=sprintf('%s %1.0d %s, ',title2, statistics(i, 1), lbl   );
                else
                    title2=sprintf('%s %1.0d %s',title2, statistics(i, 1), lbl );
                end
            end


            %2 title and labels
            if sw.alternate==0
                if isempty(mp)
                    title(axes1,{tit ...
                        title2 ...
                        });
                else
                    title(axes1,{tit ...
                        title2 ...
                        sprintf('simultanious move, n=%d, sigma=%g, eta=%g, k1=%g, k2=%g, beta=%g, tr=%g, B=(%g;%g)',mp.nC,mp.sigma, mp.eta, mp.k1, mp.k2, mp.df,mp.c_tr,mp.beta_a,mp.beta_b) ...
                        });
                end
            else
                if isempty(mp)
                    title(axes1,{tit ...
                        title2 ...
                        });
                else
                    title(axes1,{tit ...
                        title2 ...
                        sprintf('alternating move, n=%d, sigma=%g, eta=%g, k1=%g, k2=%g, beta=%g, tr=%g, B=(%g;%g), tpm11=%g tpm22=%g',mp.nC,mp.sigma, mp.eta, mp.k1, mp.k2, mp.df,mp.c_tr,mp.beta_a,mp.beta_b, mp.tpm(1,1),mp.tpm(2,2)) ...
                        });
                end
            end
        end
        
        function [statistics labels]=EqbstrPlot(eqbstr,vars,monopoly,mp,sw,tit,varargin)
            % Plots the scatterplot of the equilibrium realizations in eqbstrings
            % Assumes: eqbstring(1) - lexicographical index
            %          eqbstring(2) - count of repeated equilibria
            %          eqbstring(3) - value for x axis
            %          eqbstring(4) - value for y axis
            %          eqbstring(5,6,..) - additional statistics
            % Inputs: eqbstring - as outputed by solver
            %         vars = [i j] the columns of eqbstring to be used for size and color
            %                      if only size is given, b/w coloring is produced
            %         sw - parameters as input for solver
            %         mp - same
            %         tit - title for the graph
            %    optional handle for the figure where graphs should be
            % Output: statistics on the found equilibria
            %         cell array of some labels for the statistics
            
            %do statistics first
            labels={'total number of equilibria';'number of distinct pay-offs';'pure strategy equilibria';'symmetric equilibria';'leapfrogging equilibria';'efficiency'};
            title3={'','number of repetitions','value of firm 1','value of firm 2',labels{3:end}};
            statistics = [sum(eqbstr(:,2)) size(unique(eqbstr(:,3:4),'rows'),1) eqbstr(:,2)'*eqbstr(:,5:end)];
            %do the graph
            
            %1 prepare data and plot
            if sw.alternate==1
                %add symmetric points
                eqbstr=[eqbstr;eqbstr(:,[1 2 4 3 5:end])];
                statistics=statistics*2;
            end
            
            %make better look if number of point is not too large
            if size(eqbstr,1)<100000 && ~isempty(vars)
                %reduce some repetitions
                if (numel(vars)>1 && vars(1)==2 && vars(2)==2)
                    %if size and color are nr of repetitions
                    eqbstr=groupdata(eqbstr,[3 4],2);
                elseif (numel(vars)>1 && vars(1)==2 && vars(2)~=2)
                    %if size is nr of repetitions
                    eqbstr=groupdata(eqbstr,[3 4 vars(2)],2);
                elseif (numel(vars)==1 && vars(1)==2)
                    %if size is nr of repetitions
                    eqbstr=groupdata(eqbstr,[3 4],2);
                end
                %sort descending by size and so that bubbles are overlaping from top to bottom
                eqbstr=sortrows(eqbstr,-vars(1));
            end
            
            %color/BW
            if isempty(vars) || (numel(vars)==1 && vars(1)<=size(eqbstr,2))
                %BW case
                if ~isempty(vars)
                    title3str=sprintf('Size: %s',title3{vars(1)});
                else
                    title3str='';
                end
                if ~isempty(vars)
                    mn=min(eqbstr(:,vars(1)));
                    mx=max(eqbstr(:,vars(1)));
                    if mn==mx || isempty(vars)
                        warning 'Chosen data column does not have any variation, the size of dots is set to default'
                        if size(eqbstr,1)<100000
                            sizes=ones(size(eqbstr,1),1)*10;
                        else
                            sizes=ones(size(eqbstr,1),1)*2;
                        end
                    else
                        sizes=(eqbstr(:,vars(1))-mn)*790/(mx-mn)+8;%rescale to 10 to 700
                    end
                else
                    sizes=ones(size(eqbstr,1),1)*3;
                end
                
                %plot
                if nargin>6 && ishandle(varargin{1})
                    figure1=varargin{1};
                else
                    figure1 = figure('Color',[1 1 1],'NextPlot','new');
                end
                axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',16,'FontName','Times New Roman');
                hold(axes1,'all');
                scatter(eqbstr(:,3),eqbstr(:,4),sizes,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],'Parent',axes1);
                if ~isempty(monopoly)
                    epsi=monopoly*0.003;%gap in the monopoly line
                    line([epsi monopoly-epsi],[monopoly-epsi epsi],'Parent',axes1,'LineWidth',1,'LineStyle','--','Color',[0.25 0.25 0.25]);
                    set(axes1,'XLim',[0 monopoly]);
                    set(axes1,'YLim',[0 monopoly]);
                    set(axes1,'XTick',[get(axes1,'XTick') monopoly]);
                    set(axes1,'YTick',[get(axes1,'YTick') monopoly]);
                    a=get(axes1,'XTickLabel');
                    a(end-1,:)=' ';
                    set(axes1,'XTickLabel',a);
                else
                    maxx=max(get(axes1,'XLim'));
                    maxy=max(get(axes1,'YLim'));
                    set(axes1,'XLim',[0 min([maxx maxy])]);
                    set(axes1,'YLim',[0 min([maxx maxy])]);
                end
                
            elseif numel(vars)==2 && vars(1)<=size(eqbstr,2) && vars(2)<=size(eqbstr,2)
                %color case
                title3str=[sprintf('Size: %s',title3{vars(1)}) ' ' sprintf('Color: %s',title3{vars(2)})];
                mn=min(eqbstr(:,vars(1)));
                mx=max(eqbstr(:,vars(1)));
                if mn==mx || isempty(vars)
                    warning 'First chosen column does not have any variation, the size of dots is set to default'
                    if size(eqbstr,1)<100000
                        sizes=ones(size(eqbstr,1),1)*25;
                    else
                        sizes=ones(size(eqbstr,1),1)*5;
                    end
                else
                    sizes=(eqbstr(:,vars(1))-mn)*675/(mx-mn)+25;%rescale to 25 to 700
                end
                colors=eqbstr(:,vars(2));
                
                %plot
                if nargin>6 && ishandle(varargin{1})
                    figure1=varargin{1};
                else
                    figure1 = figure('Color',[1 1 1],'NextPlot','new');
                end
                axes1 = axes('Parent',figure1,'FontSize',14,'FontName','Times New Roman',...
                    'YGrid','on',...
                    'YColor',[0.25 0.25 0.25],...
                    'XGrid','on',...
                    'XColor',[0.25 0.25 0.25],...
                    'PlotBoxAspectRatio',[1 1 1]);
                if nargin>7 && isnumeric(varargin{2})
                    set(axes1,'CLim',[varargin{2}(1) varargin{2}(2)]);
                end
                colormap('jet');
                box(axes1,'on');
                hold(axes1,'all');
                scatter(eqbstr(:,3),eqbstr(:,4),sizes,colors,'MarkerFaceColor','flat','MarkerEdgeColor',[0 0 0],'Parent',axes1);
                colorbar('peer',axes1);
                if ~isempty(monopoly)
                    epsi=0;%gap in the monopoly line
                    line([epsi monopoly-epsi],[monopoly-epsi epsi],'Parent',axes1,'LineWidth',1,'LineStyle','-','Color',[0.25 0.25 0.25]);
                    set(axes1,'XLim',[0 monopoly]);
                    ticks=cellstr(get(axes1,'XTickLabel'));
                    ticks{end}=' ';
                    ticks{end+1}=sprintf('%4.3f',monopoly);
                    set(axes1,'XTick',[get(axes1,'XTick') monopoly]);
                    set(axes1,'XTickLabel',ticks);
                    %copy this to Y
                    set(axes1,'YLim',get(axes1,'XLim'));
                    set(axes1,'YTick',get(axes1,'XTick'));
                    set(axes1,'YTickLabel',get(axes1,'XTickLabel'));
                    
                end
                maxx=max(get(axes1,'XLim'));
                maxy=max(get(axes1,'YLim'));
                set(axes1,'XLim',[0 min(maxx,maxy)]);
                set(axes1,'YLim',[0 min(maxx,maxy)]);
            else
                error 'Unknown value of second argument (vars), need column numbers in eqbstring to be used as size and color in scatter'
            end
            
            %2 title and labels
            if sw.alternate==0
                if isempty(mp)
                    title(axes1,{tit ...
                        sprintf('%1.0d equilibria, %1.0d distinct pay-off points',statistics(1),statistics(2)) ...
                        title3str ...
                        });
                else
                    title(axes1,{tit ...
                        sprintf('%1.0d equilibria, %1.0d distinct pay-off points',statistics(1),statistics(2)) ...
                        sprintf('simultanious move, n=%d, sigma=%g, eta=%g, k1=%g, k2=%g, beta=%g, tr=%g, B=(%g;%g)',mp.nC,mp.sigma, mp.eta, mp.k1, mp.k2, mp.df,mp.c_tr,mp.beta_a,mp.beta_b) ...
                        title3str ...
                        });
                end
            else
                if isempty(mp)
                    title(axes1,{tit ...
                        sprintf('%1.0d equilibria, %1.0d distinct pay-off points',statistics(1),statistics(2)) ...
                        title3str ...
                        });
                else
                    title(axes1,{tit ...
                        sprintf('%1.0d equilibria, %1.0d distinct pay-off points',statistics(1),statistics(2)) ...
                        sprintf('alternating move, n=%d, sigma=%g, eta=%g, k1=%g, k2=%g, beta=%g, tr=%g, B=(%g;%g), tpm11=%g tpm22=%g',mp.nC,mp.sigma, mp.eta, mp.k1, mp.k2, mp.df,mp.c_tr,mp.beta_a,mp.beta_b, mp.tpm(1,1),mp.tpm(2,2)) ...
                        title3str ...
                        });
                end
            end
        end
    end
end

function res = groupdata (data, groupby, reduce)
% This function reduces the data in rows by grouping over given colums and perfoming operation on the given columns.
% Inputs: data - matrix with dara in rows, columns are variables, rows are observations
%         groupby - vector of indeces of variables to be used for grouping
%         reduce - vector of indeces of variables to be reduced using operation within each group
%         operation - pointer to function operating on columns (like sum or max)
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
end %function

