
%*************************************************************
%% ccc game
%*************************************************************
clc
clear 
!rm symbolic_sm.txt;
diary symbolic_sm.txt;
for FIRM=1:2;         % switch to sholutior for firm1 (FIRM=1) and firm2 (FIRM=2)
SIMPLIFIED=1;   % switch to show simplified version of solution (SIMPLIFIED=1)
if FIRM==1;
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('COST RECURSIONS FOR FIRM 1\n')
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('\n');

    fprintf('-----------------------------------------------\n')
    fprintf('Solving the ccc game\n')
    fprintf('-----------------------------------------------\n')
    fprintf('AN1=EPN1+beta*pc*(P2*hC1+(1-P2)*hC1)\n');
    fprintf('AI1=EPI1 +beta*pc*(P2*hC1+(1-P2)*hC1)\n'); 
    fprintf('B=beta*(1-pc)\n');

    syms AI1 AN1 B P1 P2 CI1 CN1

    S1_ccc = solve( ...
        CN1==   AN1 + B*P2    *(P1*CI1+(1-P1)*CN1) ...
                    + B*(1-P2)*(P1*CI1+(1-P1)*CN1), ... 
        CI1==   AI1 + B*P2    *(P1*CI1+(1-P1)*CN1) ...
                     +B*(1-P2)*(P1*CI1+(1-P1)*CN1), ...
       CN1,CI1);

    CN1_ccc=S1_ccc.CN1
    CI1_ccc=S1_ccc.CI1

    if ~SIMPLIFIED
    %*************************************************************
    %% c1cc game
    %*************************************************************
    fprintf('-----------------------------------------------\n')
    fprintf('Solving the c1cc game\n')
    fprintf('-----------------------------------------------\n')
    fprintf('AN1=EPN1+beta*pc*(P2*hC1+(1-P2)*hC1)\n');
    fprintf('AI1=EPI1 +beta*pc*(P2*hC1_ccc+(1-P2)*hC1_ccc)\n'); 
    fprintf('B=beta*(1-pc)\n');

    syms AI1 AN1 B P1 P2 CI1 CN1 ...
            P1_ccc CI1_ccc CN1_ccc;

    S1_c1cc = solve( ...
        CN1==   AN1 + B*P2    *(P1*CI1+(1-P1)*CN1) ...
                      + B*(1-P2)*(P1    *CI1      +(1-P1)     *CN1      ), ...
        CI1==   AI1 + B*P2    *(P1_ccc *CI1_ccc  +(1-P1_ccc)  *CN1_ccc ) ...
                       +B*(1-P2)*(P1_ccc*CI1_ccc+(1-P1_ccc) *CN1_ccc), ...
    CN1,CI1);

    CN1_c1cc=S1_c1cc.CN1
    CI1_c1cc=S1_c1cc.CI1

    else
    %*************************************************************
    %% c1cc game   - SIMPLIFIED VERSION
    %*************************************************************
    fprintf('-----------------------------------------------\n')
    fprintf('Solving the c1cc game   - SIMPLIFIED VERSION\n')
    fprintf('-----------------------------------------------\n')
    fprintf('B=beta*(1-pc)\n');
    fprintf('AN1=EPN1+beta*pc*(P2*hC1+(1-P2)*hC1)\n');
    fprintf('AI1=EPI1 +beta*pc*(P2*hC1_ccc+(1-P2)*hC1_ccc)\n'); 
    fprintf('   + B*P2    *(P1_ccc *CI1_ccc  +(1-P1_ccc)  *CN1_ccc ) \n'); 
    fprintf('   + B*(1-P2)*(P1_ccc*CI1_ccc+(1-P1_ccc) *CN1_ccc)  \n'); 

    syms AI1 AN1 B P1 P2 CI1 CN1 ...
            P1_ccc CI1_ccc CN1_ccc;

    S1_c1cc = solve( ...
        CN1==   AN1 + B*P2    *(P1*CI1+(1-P1)*CN1) ...
                      + B*(1-P2)*(P1    *CI1      +(1-P1)     *CN1      ), ...
        CI1==   AI1 , ...
    CN1,CI1);

    CN1_c1cc=S1_c1cc.CN1
    CI1_c1cc=S1_c1cc.CI1
    end

    if ~SIMPLIFIED
    %*************************************************************
    %% cc2c game 
    %*************************************************************
    fprintf('-----------------------------------------------\n')
    fprintf('Solving the cc2c game\n')
    fprintf('-----------------------------------------------\n')
    fprintf('AN1=EPN1+beta*pc*(P2*hC1_ccc+(1-P2)*hC1)\n');
    fprintf('AI1=EPI1 +beta*pc*(P2*hC1_ccc+(1-P2)*hC1)\n'); 
    fprintf('B=beta*(1-pc)\n');

    syms AI1 AN1 B P1 P2 CI1 CN1 ...
            P1_ccc CI1_ccc CN1_ccc;

    S1_cc2c = solve( ...
        CN1==   AN1 + B*P2    *(P1_ccc*CI1_ccc+(1-P1_ccc)*CN1_ccc) ...
                    + B*(1-P2)*(P1*CI1+(1-P1)*CN1), ...
        CI1==   AI1 + B*P2*(P1_ccc *CI1_ccc+(1-P1_ccc)*CN1_ccc) ...
                    + B*(1-P2)*(P1*CI1+(1-P1)*CN1), ...
    CN1,CI1);

    CN1_cc2c=S1_cc2c.CN1
    CI1_cc2c=S1_cc2c.CI1

    else
    %*************************************************************
    %% cc2c game -  SIMPLIFIED VERSION 
    %*************************************************************
    fprintf('-----------------------------------------------\n')
    fprintf('Solving the cc2c game -  SIMPLIFIED VERSION \n')
    fprintf('-----------------------------------------------\n')
    fprintf('AN1=EPN1+beta*pc*(P2*hC1_ccc+(1-P2)*hC1)\n');
    fprintf('   + B*P2    *(P1_ccc*CI1_ccc+(1-P1_ccc)*CN1_ccc)\n');
    fprintf('AI1=EPI1 +beta*pc*(P2*hC1_ccc+(1-P2)*hC1)\n'); 
    fprintf('  + B*P2*(P1_ccc *CI1_ccc+(1-P1_ccc)*CN1_ccc)\n');
    fprintf('B=beta*(1-pc)\n');

    syms AI1 AN1 B P1 P2 CI1 CN1 ...
            P1_ccc CI1_ccc CN1_ccc;

    S1_cc2c = solve( ...
        CN1==   AN1 + B*(1-P2)*(P1*CI1+(1-P1)*CN1), ...
        CI1==   AI1 + B*(1-P2)*(P1*CI1+(1-P1)*CN1), ...
    CN1,CI1);

    CN1_cc2c=S1_cc2c.CN1
    CI1_cc2c=S1_cc2c.CI1

    end

    if ~SIMPLIFIED
    %*************************************************************
    %% c1c2c game
    %*************************************************************
    fprintf('-----------------------------------------------\n')
    fprintf('Solving the c1c2c game\n')
    fprintf('-----------------------------------------------\n')
    fprintf('AN1=EPN1 +beta*pc*(P2*hC1_c1cc+(1-P2)*hC1)\n');
    fprintf('AI1=EPI1 +beta*pc*(P2*hC1_ccc+(1-P2)*hC1_cc2c)\n'); 
    fprintf('B=beta*(1-pc)\n');

    syms AI1 AN1 B P1 P2 CI1 CN1 ...
            P1_c1cc CI1_c1cc CN1_c1cc ... 
            P1_cc2c CI1_cc2c CN1_cc2c ...
            P1_ccc CI1_ccc CN1_ccc;

    S1_c1c2c = solve( ...
        CN1==   AN1 + B*P2    *(P1_c1cc*CI1_c1cc+(1-P1_c1cc)*CN1_c1cc) ...
                      + B*(1-P2)*(P1    *CI1      +(1-P1)     *CN1      ), ...
        CI1==   AI1 + B*P2    *(P1_ccc *CI1_ccc  +(1-P1_ccc)  *CN1_ccc ) ...
                       +B*(1-P2)*(P1_cc2c*CI1_cc2c+(1-P1_cc2c) *CN1_cc2c), ...
    CN1,CI1);

    CN1_c1c2c=S1_c1c2c.CN1
    CI1_c1c2c=S1_c1c2c.CI1

    else 
        
    %*************************************************************
    %% c1c2c game    - SIMPLIFIED VERSION
    %*************************************************************
    fprintf('-----------------------------------------------\n')
    fprintf('Solving the c1c2c game - SIMPLIFIED VERSION\n')
    fprintf('-----------------------------------------------\n')
    fprintf('B=beta*(1-pc)\n');
    fprintf('AN1=EPN1+beta*pc*(P2*hC1_c1cc+(1-P2)*hC1)\n');
    fprintf('   +B*P2*(P1_c1cc*CI1_c1cc+(1-P1_c1cc)*CN1_c1cc)\n\n');
    fprintf('AI1=EPI1 +beta*pc*(P2*hC1_ccc+(1-P2)*hC1_cc2c)\n')
    fprintf('   + B*P2    *(P1_ccc *CI1_ccc  +(1-P1_ccc)  *CN1_ccc )\n');    
    fprintf('   + B*(1-P2)*(P1_cc2c*CI1_cc2c+(1-P1_cc2c) *CN1_cc2c)\n')

    syms AI1 AN1 B P1 P2 CI1 CN1 ...
            P1_c1cc CI1_c1cc CN1_c1cc ... 
            P1_cc2c CI1_cc2c CN1_cc2c ...
            P1_ccc CI1_ccc CN1_ccc;

    S1_c1c2c = solve( ...
        CN1==   AN1 + B*(1-P2)*(P1    *CI1      +(1-P1)     *CN1      ), ...
        CI1==   AI1, ... 
    CN1,CI1);

    CN1_c1c2c=S1_c1c2c.CN1
    CI1_c1c2c=S1_c1c2c.CI1

    end



else  % FIRM 2
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('COST RECURSIONS FOR FIRM 2\n')
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('\n');


    %*************************************************************
    %% ccc game
    %*************************************************************
    fprintf('-----------------------------------------------\n')
    fprintf('Solving the c1c2c game - firm 2 \n')
    fprintf('-----------------------------------------------\n')
    fprintf('B=beta*(1-pc)\n');
    fprintf('AN2=EPN2+ beta*pc*(P2*hC2+(1-P2)*hC2)\n');
    fprintf('AI2=EPI2 +beta*pc*(P2*hC2+(1-P2)*hC2)\n'); 

    syms AI2 AN2 B P1 P2 CI2 CN2;
    S2_ccc = solve( ...
        CN2==   AN2 + B*P1    *(P2*CI2+(1-P2)*CN2) ...
                    + B*(1-P1)*(P2*CI2+(1-P2)*CN2), ...
        CI2==   AI2 + B*P1    *(P2*CI2+(1-P2)*CN2) ...
                    + B*(1-P1)*(P2*CI2+(1-P2)*CN2), ...
    CN2,CI2);

    CN2_ccc=S2_ccc.CN2
    CI2_ccc=S2_ccc.CI2


    %*************************************************************
    %% c1cc game
    %*************************************************************
    fprintf('-----------------------------------------------\n')
    fprintf('Solving the c1cc game - firm 2 \n')
    fprintf('-----------------------------------------------\n')
    fprintf('B=beta*(1-pc)\n');
    fprintf('AN2=EPN2+ beta*pc*(P2*hC2_ccc+(1-P2)*hC2)\n');
    fprintf('+ B*P1*(P2_ccc*CI2_ccc+(1-P2_ccc)*CN2_ccc)\n');
    fprintf('AI2=EPI2 +beta*pc*(P2*hC2_ccc+(1-P2)*hC2)\n'); 
    fprintf('   +B*P1*(P2_ccc*CI2_ccc+(1-P2_ccc)*CN2_ccc )\n');

    syms AI2 AN2 B P1 P2 CI2 CN2 ...
            P2_ccc  CI2_ccc  CN2_ccc;

    S2_c1cc = solve( ...
        CN2==   AN2 + B*(1-P1)*(P2    *CI2      +(1-P2)     *CN2      ), ...
        CI2==   AI2 + B*(1-P1)*(P2*CI2+(1-P2)*CN2), ...
    CN2,CI2);

    CN2_c1cc=S2_c1cc.CN2
    CI2_c1cc=S2_c1cc.CI2


    %*************************************************************
    %% cc2c game
    %*************************************************************
    fprintf('-----------------------------------------------\n')
    fprintf('Solving the cc2c game - firm 2 \n')
    fprintf('-----------------------------------------------\n')
    fprintf('B=beta*(1-pc)\n');
    fprintf('AN2=EPN2+ beta*pc*(P2*hC2_cc2c+(1-P2)*hC2)\n');
    fprintf('AI2=EPI2 +beta*pc*(P2*hC2_ccc+(1-P2)*hC2_ccc)\n'); 
    fprintf('  + B*P1    *(P2_ccc*CI2_ccc  +(1-P2_ccc) *CN2_ccc ) \n'); 
    fprintf('  + B*(1-P1)*(P2_ccc*CI2_ccc+(1-P2_ccc)*CN2_ccc)\n'); 

    syms AI2 AN2 B P1 P2 CI2 CN2 ...
            P2_ccc  CI2_ccc  CN2_ccc;

    S2_cc2c = solve( ...
        CN2==   AN2 + B*P1    *(P2*CI2+(1-P2)*CN2) ...
                    + B*(1-P1)*(P2*CI2+(1-P2)*CN2), ...
        CI2==   AI2 , ...
    CN2,CI2);

    CN2_cc2c=S2_cc2c.CN2
    CI2_cc2c=S2_cc2c.CI2


    if SIMPLIFIED
        %*************************************************************
        %% c1c2c game - Simplified
        %*************************************************************
        fprintf('-----------------------------------------------\n')
        fprintf('Solving the c1c2c game - firm 2 \n')
        fprintf('-----------------------------------------------\n')
        fprintf('B=beta*(1-pc)\n');
        fprintf('AN2=EPN2+ beta*pc*(P2*hC2_cc2c+(1-P2)*hC2)\n');
        fprintf('+ B*P1    *(P2_cc2c*CI2_cc2c+(1-P2_cc2c)*CN2_cc2c)\n'); 
        fprintf('AI2=EPI2 +beta*pc*(P2*hC2_ccc+(1-P2)*hC2_c1cc)\n'); 
        fprintf('+ B*P1    *(P2_ccc*CI2_ccc  +(1-P2_ccc) *CN2_ccc )\n');
        fprintf('    + B*(1-P1)*(P2_c1cc*CI2_c1cc+(1-P2_c1cc)*CN2_c1cc)\n');

        syms AI2 AN2 B P1 P2 CI2 CN2

        S2_c1c2c = solve( ...
            CN2==   AN2 + B*(1-P1)*(P2    *CI2      +(1-P2)     *CN2      ), ...
            CI2==   AI2,  ...
        CN2,CI2);

        CN2_c1c2c=S2_c1c2c.CN2
        CI2_c1c2c=S2_c1c2c.CI2
    else
        %*************************************************************
        %% c1c2c game - NOT simplified
        %*************************************************************
        fprintf('-----------------------------------------------\n')
        fprintf('Solving the c1c2c game - firm 2 - NOT SIMPLIFIED \n')
        fprintf('-----------------------------------------------\n')
        fprintf('B=beta*(1-pc)\n');
        fprintf('AN2=EPN2+ beta*pc*(P2*hC2_cc2c+(1-P2)*hC2)\n');
        fprintf('AI2=EPI2 +beta*pc*(P2*hC2_ccc+(1-P2)*hC2_c1cc)\n'); 

        syms AI2 AN2 B P1 P2 CI2 CN2 ...
                P2_c1cc CI2_c1cc CN2_c1cc ... 
                P2_cc2c CI2_cc2c CN2_cc2c ...
                P2_ccc  CI2_ccc  CN2_ccc;

        S2_c1c2c = solve( ...
            CN2==   AN2 + B*P1    *(P2_cc2c*CI2_cc2c+(1-P2_cc2c)*CN2_cc2c) ...
                        + B*(1-P1)*(P2    *CI2      +(1-P2)     *CN2      ), ...
            CI2==   AI2 + B*P1    *(P2_ccc*CI2_ccc  +(1-P2_ccc) *CN2_ccc ) ...
                        + B*(1-P1)*(P2_c1cc*CI2_c1cc+(1-P2_c1cc)*CN2_c1cc), ...
        CN2,CI2);

        CN2_c1c2c=S2_c1c2c.CN2
        CI2_c1c2c=S2_c1c2c.CI2
    end  % End SIMPLIFIED
end  % End FIRM=1,2
end % of loop over firms
diary off;




