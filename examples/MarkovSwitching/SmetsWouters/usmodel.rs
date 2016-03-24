endogenous ewma epinfma  zcapf rkf kf pkf cf invef yf labf wf rrf mc zcap rk k pk c inve y lab pinf w r a  b g
qs  ms  spinf sw kpf kp 

endogenous labobs robs pinfobs dy dc dinve dw
 
exogenous ea eb eg  eqs  em  epinf ew  
 
parameters curvw cgy curvp constelab constepinf constebeta cmaw cmap calfa 
czcap csadjcost ctou csigma chabb cfc 
cindw cprobw cindp cprobp csigl clandaw 
crhoa crhoas crhob crhog crhols crhoqs crhoms crhopinf crhow  
ctrend cg, sig_a, sig_b,sig_w,sig_pinf,sig_m,sig_qs,sig_g

parameters(coef,2)  crpi crdy cry crr

parameters coef_tp_1_2, coef_tp_2_1

model 

		# cpie=1+constepinf/100;
		# cgamma=1+ctrend/100 ;
		# cbeta=1/(1+constebeta/100);
		# clandap=cfc;
		# cbetabar=cbeta*cgamma^(-csigma);
		# cr=cpie/(cbeta*cgamma^(-csigma));
		# crk=(cbeta^(-1))*(cgamma^csigma) - (1-ctou);
		# cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));
		# cikbar=(1-(1-ctou)/cgamma);
		# cik=(1-(1-ctou)/cgamma)*cgamma;
		# clk=((1-calfa)/calfa)*(crk/cw);
		# cky=cfc*(clk)^(calfa-1);
		# ciy=cik*cky;
		# ccy=1-cg-cik*cky;
		# crkky=crk*cky;
		# cwhlc=(1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
		# cwly=1-crk*cky;
		# conster=(cr-1)*100;

% flexible economy

	      0*(1-calfa)*a + 1*a =  calfa*rkf+(1-calfa)*(wf)  ;
	      zcapf =  (1/(czcap/(1-czcap)))* rkf  ;
	      rkf =  (wf)+labf-kf ;
	      kf =  kpf(-1)+zcapf ;
	      invef = (1/(1+cbetabar*cgamma))* (  invef(-1) + cbetabar*cgamma*invef(1)+(1/(cgamma^2*csadjcost))*pkf ) +qs ;
          pkf = -rrf-0*b+(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b +(crk/(crk+(1-ctou)))*rkf(1) +  ((1-ctou)/(crk+(1-ctou)))*pkf(1) ;
	      cf = (chabb/cgamma)/(1+chabb/cgamma)*cf(-1) + (1/(1+chabb/cgamma))*cf(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(labf-labf(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(rrf+0*b) + b ;
	      yf = ccy*cf+ciy*invef+g  +  crkky*zcapf ;
	      yf = cfc*( calfa*kf+(1-calfa)*labf +a );
	      wf = csigl*labf 	+(1/(1-chabb/cgamma))*cf - (chabb/cgamma)/(1-chabb/cgamma)*cf(-1) ;
	      kpf =  (1-cikbar)*kpf(-1)+(cikbar)*invef + (cikbar)*(cgamma^2*csadjcost)*qs ;

% sticky price - wage economy

	      mc =  calfa*rk+(1-calfa)*(w) - 1*a - 0*(1-calfa)*a ;
	      zcap =  (1/(czcap/(1-czcap)))* rk ;
	      rk =  w+lab-k ;
	      k =  kp(-1)+zcap ;
	      inve = (1/(1+cbetabar*cgamma))* (  inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*pk ) +qs ;
          pk = -r+pinf(1)-0*b +(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b + (crk/(crk+(1-ctou)))*rk(1) +  ((1-ctou)/(crk+(1-ctou)))*pk(1) ;
	      c = (chabb/cgamma)/(1+chabb/cgamma)*c(-1) + (1/(1+chabb/cgamma))*c(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) + 0*b) +b ;
	      y = ccy*c+ciy*inve+g  +  1*crkky*zcap ;
	      y = cfc*( calfa*k+(1-calfa)*lab +a );
	      pinf =  (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*pinf(1) +cindp*pinf(-1) 
               +((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc)  )  + spinf ; 
	      w =  (1/(1+cbetabar*cgamma))*w(-1)
               +(cbetabar*cgamma/(1+cbetabar*cgamma))*w(1)
               +(cindw/(1+cbetabar*cgamma))*pinf(-1)
               -(1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf
               +(cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf(1)
               +(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*
               (csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1) -w) 
               + 1*sw ;

		  % Monetary policy reaction function
		  r =  crpi*(1-crr)*pinf +cry*(1-crr)*(y-yf)+crdy*(y-yf-y(-1)+yf(-1))+crr*r(-1)+ms  ;

		  a = crhoa*a(-1)  + sig_a*ea;	%
	      b = crhob*b(-1) + sig_b*eb;
	      g = crhog*(g(-1)) + sig_g*eg + cgy*sig_a*ea;
	      qs = crhoqs*qs(-1) + sig_qs*eqs;
	      ms = crhoms*ms(-1) + sig_m*em;
	      spinf = crhopinf*spinf(-1) + epinfma - cmap*epinfma(-1);
	          epinfma=sig_pinf*epinf;
	      sw = crhow*sw(-1) + ewma - cmaw*ewma(-1) ;
	      ewma=sig_w*ew; 
	      kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs ;

% measurment equations

		dy=y-y(-1)+ctrend;
		dc=c-c(-1)+ctrend;
		dinve=inve-inve(-1)+ctrend;
		dw=w-w(-1)+ctrend;
		pinfobs = 1*(pinf) + constepinf;
		robs =    1*(r) + conster;
		labobs = lab + constelab;


parameterization
	coef_tp_1_2 ,   0.1000;
	coef_tp_2_1 ,	0.7000;
	crpi(coef,1),   1.4880;
	crr(coef,1) ,   0.8762;
	cry(coef,1) ,   0.0593;
	crdy(coef,1),   0.2347;
	% zero-lower-bound-financial-crisis time
	crpi(coef,2),   0.0000;
	crr(coef,2) ,   0.0000;
	cry(coef,2) ,   0.0000;
	crdy(coef,2),   0.0000;
	
	crhoas,1; 	 % this parameter does not enter the model
	crhols,    0.9928;  % this parameter does not enter the model
	ctou,.025;
	clandaw,1.5;
	cg,0.18;
	curvp,10;
	curvw,10;
	
	sig_a,   0.4618  ;
	sig_b,   0.18185 ;
	sig_g,   0.6090  ;
	sig_qs,  0.46017 ;
	sig_m,   0.2397  ;
	sig_pinf,0.1455  ;
	sig_w,   0.2089  ;
	
	calfa,   .24     ;
	csigma,   1.5    ;
	cfc,      1.5    ;
	cgy,      0.51   ;
	
	csadjcost, 6.0144;
	chabb,     0.6361;
	cprobw,    0.8087;
	csigl,     1.9423;
	cprobp,    0.6   ;
	cindw,     0.3243;
	cindp,     0.47  ;
	czcap,     0.2696;
	crhoa,     0.9977;
	crhob,     0.5799;
	crhog,     0.9957;
	crhoqs,    0.7165;
	crhoms,    .3    ;
	crhopinf,  0.8   ;
	crhow,     0.8   ;
	cmap ,     0.7   ;
	cmaw  ,    0.7   ;
	constebeta,0.7420;
	
	ctrend,    0.3982;
	constepinf, .7   ;
	constelab, 1.2918;
