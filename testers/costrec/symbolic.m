%% DEFINE vars
clc
clear
fprintf('DEFINE\n');
fprintf('-----------------------------------------------\n')
fprintf('A1= f11*HC1 + f21*HC2 \n');
fprintf('A2= f12*HC1 + f22*HC2 \n');
fprintf('A1_ccc= f11*HC1_ccc + f21*HC2_ccc \n');
fprintf('A2_ccc= f12*HC1_ccc + f22*HC2_ccc \n');
fprintf('B1= bet*(1-pc)*f11\n');
fprintf('B2= bet*(1-pc)*f12\n');
fprintf('C1= bet*(1-pc)*f21\n');
fprintf('C2= bet*(1-pc)*f22\n');
fprintf('\n')

%% CCC
clear 
fprintf('Solving the CCC game\n')
fprintf('-----------------------------------------------\n')
syms EP11 EP10 EPx1 CN1 CI1 CI2 CN2 A1 A2 B1 B2 C1 C2  P1 P2

S_ccc = solve(...
    CN1 == EP10 + A1 + B1*(P1*CI1 + (1-P1)*CN1) + C1*(P2*CI2 + (1-P2)*CN2), ...
    CI1 == EP11 + A1 + B1*(P1*CI1 + (1-P1)*CN1) + C1*(P2*CI2 + (1-P2)*CN2), ...
    CN2 == EPx1 + A2 + B2*(P1*CI1 + (1-P1)*CN1) + C2*(P2*CI2 + (1-P2)*CN2), ...
    CI2 == EPx1 + A2 + B2*(P1*CI1 + (1-P1)*CN1) + C2*(P2*CI2 + (1-P2)*CN2), ...
    CN1,CI1, CN2, CI2);
 
CI2_ccc=S_ccc.CI2


S_ccc = solve(...
    CN1 == EP10 + A1 + B1*(P1*CI1 + (1-P1)*CN1) + C1*(P2*CI2 + (1-P2)*CN2), ...
    CI1 == EP11 + A1 + B1*(P1*CI1 + (1-P1)*CN1) + C1*(P2*CI2 + (1-P2)*CN2), ...
    CN2 == EPx1 + A2 + B2*(P1*CI1 + (1-P1)*CN1) + C2*(P2*CI2 + (1-P2)*CN2), ...
    CN1,CI1, CN2);

CN1_ccc=S_ccc.CN1

S_ccc = solve(...
    CN2 == CI2, ...
    CI1 == EP11 + A1 + B1*(P1*CI1 + (1-P1)*CN1) + C1*(P2*CI2 + (1-P2)*CN2), ...
    CI1, CN2);

CN2_ccc=S_ccc.CN2
CI1_ccc=S_ccc.CI1

%    CN2 == EPx + A2 + B2*(P1*CI1 + (1-P1)*CN1) + C2*(P2*CI2 + (1-P2)*CN2), ...



fprintf('\n')

%% solving the C1CC game
clear 
fprintf('Solving the C1CC game\n')
fprintf('-----------------------------------------------\n')

syms EP11 EP10 EPx1 CN1_ccc CI1_ccc CI2_ccc CN2_ccc A1 A2 A1_ccc A2_ccc B1 B2 C1 C2  P1 P2

CI1_c1cc= EP11 + A1_ccc + B1*(P1*CI1_ccc + (1-P1)*CN1_ccc) + C1*(P2*CI2_ccc + (1-P2)*CN2_ccc)

syms EP11 EP10 EPx1 CN1 CI1 CI2 CN2 A1 A2 B1 B2 C1 C2  P1 P2

S_c1cc = solve(...
    CN1 == EP10 +A1 + B1*(P1*CI1 + (1-P1)*CN1) + C1*(P2*CI2 + (1-P2)*CN2), ...
    CN2 == EPx1 +A2 + B2*(P1*CI1 + (1-P1)*CN1) + C2*(P2*CI2 + (1-P2)*CN2), ...
    CI2 == EPx1 +A2 + B2*(P1*CI1 + (1-P1)*CN1) + C2*(P2*CI2 + (1-P2)*CN2), ...
    CN1, CN2, CI2);
 
CN1_c1cc=S_c1cc.CN1

% S_c1cc = solve(...
%     CN2 == EPx1 +A2 + B2*(P1*CI1 + (1-P1)*CN1) + C2*(P2*CI2 + (1-P2)*CN2), ...
%     CI2 == EPx1 +A2 + B2*(P1*CI1 + (1-P1)*CN1) + C2*(P2*CI2 + (1-P2)*CN2), ...
%     CN2, CI2);
% 
CN2_c1cc=S_c1cc.CN2
CI2_c1cc=S_c1cc.CI2


%% solving the CC2C game
clear 
fprintf('Solving the CC2C game\n')
fprintf('-----------------------------------------------\n')

syms EP11 EP10 EPx1 CN1_ccc CI1_ccc CI2_ccc CN2_ccc A1 A2 A1_ccc A2_ccc B1 B2 C1 C2  P1 P2

CI2_cc2c= EPx1 + A2_ccc + B2*(P1*CI1_ccc + (1-P1)*CN1_ccc) + C2*(P2*CI2_ccc + (1-P2)*CN2_ccc)
syms EP11 EP10 EPx1 CN1 CI1 CI2 CN2 A1 A2 B1 B2 C1 C2  P1 P2

S_cc2c = solve(...
    CN1 == EP10 + A1 + B1*(P1*CI1 + (1-P1)*CN1) + C1*(P2*CI2 + (1-P2)*CN2), ...
    CI1 == EP11 + A1 + B1*(P1*CI1 + (1-P1)*CN1) + C1*(P2*CI2 + (1-P2)*CN2), ...
    CN2 == EPx1  + A2 + B2*(P1*CI1 + (1-P1)*CN1) + C2*(P2*CI2 + (1-P2)*CN2), ...
    CN1,CI1, CN2);
 
CN1_cc2c=S_cc2c.CN1
CI1_cc2c=S_cc2c.CI1
CN2_cc2c=S_cc2c.CN2

%% Expected cost recursions for c1,c2,c game for firm 1 
% CN1 Expected cost of not investing, when firm 1 moves
% CN2 Expected cost of not investing, when firm 2 moves

clear
fprintf('Solving the C1C2C game\n')
fprintf('-----------------------------------------------\n')
syms A1 B1 C1 A2 B2 C2 CN1 CN2
fCN1(A1, B1, C1, CN1, CN2)  = A1 + B1*CN1 + C1*CN2;
fCN2(A2, B2, C2, CN1, CN2)  = A2 + B2*CN1 + C2*CN2;

CI1_c1c2c='closed form'
CI2_c1c2c='closed form'

Sc1c2c = solve(...
    CN1 == fCN1(A1, B1, C1, CN1, CN2), ...
    CN2 == fCN2(A2, B2, C2, CN1, CN2) , ...
    CN1,CN2);
CN1_c1c2c=Sc1c2c.CN1

Sc1c2c = solve(...
    CN2 == fCN2(A2, B2, C2, CN1, CN2) , ...
    CN2);
CN2_c1c2c=Sc1c2c
return

%% FIRM 2!!!

%% solving the C1CC game FOR FIRM2
clc
clear 
fprintf('Solving the c1cc game\n FIRM 2')
fprintf('-----------------------------------------------\n')

syms EP21 EP20 EPx2 CVN CVI CXI CXN A1 A2 B1 B2 C1 C2  P1 P2

S_c1cc = solve(...
    CVN == EP20 + A1 + B1*(P2*CVI + (1-P2)*CVN) + C1*(P1*CXI + (1-P1)*CXN), ...
    CVI == EP21 + A1 + B1*(P2*CVI + (1-P2)*CVN) + C1*(P1*CXI + (1-P1)*CXN), ...
    CXN == EPx2 + A2 + B2*(P2*CVI + (1-P2)*CVN) + C2*(P1*CXI + (1-P1)*CXN), ...
    CVN,CVI, CXN);
 
CXI='closed form'
CVN_c1cc=S_c1cc.CVN
CVI_c1cc=S_c1cc.CVI
CXN_c1cc=S_c1cc.CXN

%% solving the CC2C game - for firm 2
clear 
fprintf('Solving the cc2c game for frim 2\n')
fprintf('-----------------------------------------------\n')

syms EP21 EP20 EPx1 CVN_ccc CVI_ccc CXI_ccc CXN_ccc A1 A2 A1_ccc A2_ccc B1 B2 C1 C2  P1 P2

CVI_cc2c= EP21 + A1_ccc + B1*(P2*CVI_ccc + (1-P2)*CVN_ccc) + C1*(P1*CXI_ccc + (1-P1)*CXN_ccc)

syms EP21 EP20 EPx2 CVN CVI CXI CXN A1 A2 B1 B2 C1 C2  P1 P2

S_cc2c = solve(...
    CVN == EP20 +A1 + B1*(P2*CVI + (1-P2)*CVN) + C1*(P1*CXI + (1-P1)*CXN), ...
    CXN == EPx2 +A2 + B2*(P2*CVI + (1-P2)*CVN) + C2*(P1*CXI + (1-P1)*CXN), ...
    CXI == EPx2 +A2 + B2*(P2*CVI + (1-P2)*CVN) + C2*(P1*CXI + (1-P1)*CXN), ...
    CVN, CXN, CXI);
 
CVN_cc2c=S_cc2c.CVN
CXN_cc2c=S_cc2c.CXN
CXI_cc2c=S_cc2c.CXI

CVI_cc2c =
 
A1_ccc + EP21 - B1*(CVN_ccc*(P2 - 1) - CVI_ccc*P2) - C1*(CXN_ccc*(P1 - 1) - CXI_ccc*P1)
 
 
CVN_cc2c =
 
-(A1 + EP20 - A1*C2 + A2*C1 - C2*EP20 + C1*EPx2 + B1*CVI*P2 - B1*C2*CVI*P2 + B2*C1*CVI*P2)/(B1 + C2 - B1*C2 + B2*C1 - B1*P2 + B1*C2*P2 - B2*C1*P2 - 1)
 
 
CXN_cc2c =
 
-(A2 + EPx2 + A1*B2 - A2*B1 + B2*EP20 - B1*EPx2 - A1*B2*P2 + A2*B1*P2 + B2*CVI*P2 - B2*EP20*P2 + B1*EPx2*P2)/(B1 + C2 - B1*C2 + B2*C1 - B1*P2 + B1*C2*P2 - B2*C1*P2 - 1)
 
 
CXI_cc2c =
 
-(A2 + EPx2 + A1*B2 - A2*B1 + B2*EP20 - B1*EPx2 - A1*B2*P2 + A2*B1*P2 + B2*CVI*P2 - B2*EP20*P2 + B1*EPx2*P2)/(B1 + C2 - B1*C2 + B2*C1 - B1*P2 + B1*C2*P2 - B2*C1*P2 - 1)
 






