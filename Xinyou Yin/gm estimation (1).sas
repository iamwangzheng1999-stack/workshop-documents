TITLE 'a NRH-A model to estimate measophyll conductance';
DATA A;
INPUT O IINC CI A PhI2 ;
CARDS;
21 1000.1 43.967 0.122 0.2033
21 1000.1 75.614 2.706 0.2351
21 1000.1 142.740 10.971 0.2839
;
PROC NLIN DATA=A METHOD=GAUSS ITERATIONS=100;
PARMS GM=0.45;

LUMP = 0.4543;
RD = 0.6195;

SCO = 3.022;
GAMMAX=0.5*(O*10)/SCO;

X1 = LUMP*IINC*PHI2/4;

BBJ = (X1-RD)/GM+2*GAMMAX+CI;
SQJ = (BBJ**2-4/GM*((CI-GAMMAX)*X1-(2*GAMMAX+CI)*RD))**0.5;

MODEL A = (BBJ-SQJ)/(2/GM);

output out = b predicted = yp residual = res ;
proc corr;
var A yp;
proc print;
RUN;
