% Rudebusch-Swanson model by Andrew Binning
% "The Bond Premium in a DSGE Model with Long-Run and Nominal Risks"
% By Glenn D. Rudebusch and Eric T. Swanson. American Economic Journal:
% Macroeconomics 2012, 4(1): 105-143.

endogenous V, C, L, VK, PI, Valphaexp, ZN, MC, Y, ZD, P0, W, Int, Delta, PIAVG, G, A, DZ, PISTAR,
	PB, PBB, YTM, YTMM, TT

exogenous EPSInt, EPSA, EPSZ, EPSG, EPSPISTAR

parameters phi, chi, theta, eta, xi, alpha, LMax, DZBAR, YBAR, delta, KRAT, GRAT , 
    PIBAR, rhoinflavg, taylrho, taylpi, tayly, rhoa, rhoz, rhog, rhopistar, GSSLOAD,
	cdelta Beta_tilde
	% stand-ins for shocks' standard deviations
	%-------------------------------------------
	std_EPSInt, std_EPSA, std_EPSZ, std_EPSG, std_EPSPISTAR
	% parameters determined in steady state
	%---------------------------------------
	KBAR, IBAR chi0, Beta

model

    % Behavioural Equations
    %-----------------------
	
    V - C^(1-phi)/(1-phi) - chi0*(LMax - L)^(1-chi)/(1-chi) - Beta*VK;
	
    C^(-phi) - Beta*(exp(Int)/PI(1))*C(1)^(-phi)*DZ(1)^(-phi)*(V(1)*DZ(1)^(1-phi)/VK)^(-alpha);
	
    Valphaexp - ((V(1)/steady_state(V))*(DZ(1)/DZBAR)^(1-phi))^(1-alpha);
	
    VK - steady_state(V)*DZBAR^(1-phi)*Valphaexp^(1/(1-alpha));
	
    ZN - (1+theta)*MC*Y - xi*Beta*(C(1)/C)^(-phi)*DZ(1)^(-phi)*(V(1)*DZ(1)^(1-phi)/VK)^(-alpha)*PI(1)^((1+theta)/theta/eta)*ZN(1);
	
    ZD - Y - xi*Beta*(C(1)/C)^(-phi)*DZ(1)^(-phi)*(V(1)*DZ(1)^(1-phi)/VK)^(-alpha)*PI(1)^(1/theta)*ZD(1);
	
    P0^(1+(1+theta)/theta*(1-eta)/eta) - ZN/ZD;
	
    PI^(-1/theta) - (1-xi)*(P0*PI)^(-1/theta) - xi;
	
    MC - W/eta*Y^((1-eta)/eta)/A^(1/eta)/KBAR^((1-eta)/eta);
	
    chi0*(LMax - L)^(-chi)/C^(-phi) - W;
	
    Y - A*KBAR^(1-eta)*L^eta/Delta;
	
    Delta^(1/eta) - (1-xi)*P0^(-(1+theta)/theta/eta) - xi*PI^((1+theta)/theta/eta)*Delta(-1)^(1/eta);
	
    C - Y + G + IBAR;
	
    log(PIAVG) - rhoinflavg*log(PIAVG(-1)) - (1-rhoinflavg)*log(PI);
	
    4*Int - (1 - taylrho)*(4*log(1/Beta*DZBAR^phi) + 4*log(PIAVG) + taylpi*(4*log(PIAVG) - PISTAR) + tayly*(Y-YBAR)/YBAR ) - taylrho*4*Int(-1) - std_EPSInt*EPSInt;
	
    % Explaining the term premia calculations (next 5 equations)
    %------------------------------------------------------------
	% Rudebusch and Swanson calculate the Term Premium on a 10 year bond with the 10 year Treasury note as the bench-mark bond. To simplify calculations, 
	% they assume that households can buy and sell a long-term default free nominal consol which pays a geometrically declining coupon in every period in
	% perpetuity. The price of the risky bond is PB in period t, the price of the Treasury bond is PBB in period t.  cdelta is the rate of decay of the 
	% coupon on the consol. YTM and YTMM are the continuously compounded yields to maturity for the risky and the less risky bonds. The term premium, TT, 
	% is calculated as the difference between the yields on the risky and the less risky bonds. 
	
    PB - 1 - PB(1)*cdelta*Beta*(C(1)/C)^(-phi)*DZ(1)^(-phi)*(V(1)*DZ(1)^(1-phi)/VK)^(-alpha)/PI(1);
	
    PBB - 1 - PBB(1)*cdelta/exp(Int);
	
    YTM - log(cdelta*PB/(PB-1))*400;
	
    YTMM - log(cdelta*PBB/(PBB-1))*400;
	
    TT - 100*(YTM-YTMM);
	
    % Shock processes
    %-----------------
	
    log(A/steady_state(A)) - rhoa*log(A(-1)/steady_state(A)) - std_EPSA*EPSA;
	
    log(DZ/DZBAR) - rhoz*log(DZ(-1)/DZBAR) - std_EPSZ*EPSZ;
	
    log(G/steady_state(G)) - rhog*log(G(-1)/steady_state(G)) - std_EPSG*EPSG;
	
    PISTAR - (1-rhopistar)*PIBAR + rhopistar*PISTAR(-1) - GSSLOAD*(4*log(PIAVG) - PISTAR) - std_EPSPISTAR*EPSPISTAR;

steady_state_model
	L = 1/3;
	
	Y = 1;
	
	MC = 1/(1+theta);
	
	DZ = DZBAR;
	
	Beta = Beta_tilde*DZ^phi;
	
	PI = exp(PIBAR);
	
	Int = log((PI/Beta)*DZ^phi);
	
	Delta = 1;
	
	PIAVG = PI;
	
	KBAR = 4*YBAR*KRAT;
	
	IBAR = delta*KBAR;
	
	G = YBAR*GRAT;
	
	A =  YBAR/(KBAR^(1-eta)*L^eta);
	
	PISTAR = PIBAR;
	
	C = Y - G - IBAR;
	
	ZN = (1+theta)*MC*Y/(1-xi*Beta*DZ^(-phi)*PI^(eta*(1+theta)/theta));
	
	ZD = Y/(1-xi*Beta*DZ^(-phi)*PI^(1/theta));
	
	P0 = (ZN/ZD)^(1/(1+(1+theta)/theta*(1-eta)/eta));
	
	W = MC*eta*Y/L;
	
	chi0 = W/((LMax - L)^(-chi)/C^(-phi));
	
	xx_ssmdef_1 = C^(1-phi)/(1-phi) + chi0*(LMax - L)^(1-chi)/(1-chi);
	
	V = xx_ssmdef_1/(1 - (Beta*DZBAR^(1-phi)));
	
	VK = V*DZBAR^(1-phi);
	
	Valphaexp = (V/V)^(1-alpha);
	
	PB = 1/(1 - cdelta/exp(Int));
	
	PBB = 1/(1 - cdelta/exp(Int));
	
	YTM = log(cdelta*PB/(PB-1))*400;
	
	YTMM = log(cdelta*PBB/(PBB-1))*400;


parameterization
	phi , 2;
	chi , 3;
	theta , 0.2;
	eta , 2/3;
	xi , 0.75;
	alpha , -148.3;
	LMax , 1;
	DZBAR , 1.0025;
	YBAR , 1;
	delta , 0.02;
	KRAT , 2.5;
	GRAT , 0.17;
	PIBAR , 0;
	rhoinflavg , 0;
	taylrho , 0.73;
	taylpi , 0.53;
	tayly , 0.93;
	rhoa , 0.95;
	rhoz , 0;
	rhog , 0.95;
	rhopistar , 0;
	GSSLOAD , 0;
	Beta_tilde , 0.99;
	cdelta , 0.9848;
	std_EPSInt   , 0.01;
	std_EPSA     , 0.01;
	std_EPSZ        , 0.01;
	std_EPSG       , 0.01;
	std_EPSPISTAR , 0.01;
