%*************************************************************
%% ccc game
%*************************************************************
clc
clear 
!rm symbolic_am_solution.txt;
diary symbolic_am_solution.txt;
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf('Solving the alternativ move game\n')
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf('\n');

fprintf('-----------------------------------------------\n')
fprintf('Solving the ccc game\n')
fprintf('-----------------------------------------------\n')
fprintf('AI1v = r1-Kc + bet*pc*H11+bet*(1-pc)*tpm11*logsumK\n');
fprintf('AI1x=  r1    + bet*pc*H12+bet*(1-pc)*tpm12*logsumK\n'); 
fprintf('Bij=beta*(1-pc)*tpmij\n');


syms AI1x AI1v B11 B21 B12 B22 v11 x11 


S1_ccc = solve( ...
	v11==AI1v+(B11*v11+B21*x11), ... 
    x11==AI1x+(B12*v11+B22*x11), ... 
    x11, v11)

v11_ccc=S1_ccc.v11
x11_ccc=S1_ccc.x11

diary off;


return

S1_ccc = solve( ...
	CN1==   AN1 + B*P2    *(P1*CI1+(1-P1)*CN1) ...
	+ B*(1-P2)*(P1*CI1+(1-P1)*CN1), ... 
	CI1==   AI1 + B*P2    *(P1*CI1+(1-P1)*CN1) ...
	+B*(1-P2)*(P1*CI1+(1-P1)*CN1), ...
	CN1,CI1);

CN1_ccc=S1_ccc.CN1
CI1_ccc=S1_ccc.CI1

    %*************************************************************
    %% c1cc game  
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

    %*************************************************************
    %% cc2c game 
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









