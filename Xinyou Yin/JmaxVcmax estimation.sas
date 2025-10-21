TITLE 'Combined FvCB-photosynthesis and mesophyll conductance model';
DATA A;
   INPUT O IINC CI A ;
   CARDS;
21 1000.1 43.967 0.122
21 1000.1 75.614 2.706
21 1000.1 142.74 10.971
21 1000 292.806 20.242
21 1000.2 516.207 25.201
21 1000 851.093 26.67
21 1000 1335.018 27.459			
2 199.9 903.592 11.573
2 150 911.601 8.874
2 120.2 913.756 7.072
2 90 920.126 5.122
2 59.9 928.752 3.316
2 40.3 936.671 2.189
;
PROC NLIN DATA=A METHOD=GAUSS ITERATIONS=1000;
PARMS VCMAX=50 JMAX=200 THETA=0.5;

Rd   = 0.6183;
LUMP = 0.4542;
GM   = 0.1945;
F2LL = 0.6101;

*** Rubisco-parameters a priori assigned (based on Cousins et al. 2010);
      KMC = 291;  
      KMO = 194;
      SCO = 3.002;
      GAMMAX = 0.5*(O*10)/SCO;

*** Rubisco-limited part;
 	   X1R   = VCMAX;
       X2R   = KMC*(1+(O*10)/KMO);
       
       BBR   = (X1R-RD)+GM*(CI+X2R);
       CCR   = GM*(X1R*(CI-GAMMAX)-RD*(CI+X2R));
       
       AR     = (BBR-SQRT(BBR**2-4.*CCR))/2;

*** Electron transport limited part;
       K2LL  = LUMP*F2LL;
       BB    = K2LL*IINC + JMAX;
	   J     = (BB-SQRT(BB**2-4*THETA*JMAX*K2LL*IINC))/(2*THETA);

       X1J   = J/4;
       X2J   = 2.*GAMMAX;
       
       BBJ   = (X1J-RD)+GM*(CI+X2J);
       CCJ   = GM*(X1J*(CI-GAMMAX)-RD*(CI+X2J));

       AJ     = (BBJ-SQRT(BBJ**2-4.*CCJ))/2;

MODEL  A     = MIN(AR,AJ);

output out = b predicted = yp residual = res ;
proc corr;
 var A yp;
proc print;
RUN;
