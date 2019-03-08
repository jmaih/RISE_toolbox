endogenous ewma epinfma  zcapf rkf kf pkf cf invef yf labf wf rrf mc zcap rk k pk c inve y lab pinf w r a
b g qs  ms  spinf sw kpf kp

endogenous labobs robs pinfobs dy dc dinve dw

endogenous dr pinf_target
 
exogenous ea eb eg  eqs  em  epinf ew   
 
parameters curvw cgy curvp constelab constepinf constebeta cmaw cmap calfa 
czcap csadjcost ctou csigma chabb cfc 
cindw cprobw cindp cprobp csigl clandaw  crpi crdy cry crr 
crhoa crhoas crhob crhog crhols crhoqs crhoms crhopinf crhow  
ctrend cg, sig_a, sig_b,sig_w,sig_pinf,sig_m,sig_qs,sig_g

% loss function weights
parameters wy, wr

% we specify the probability of commitment (if it is to vary)
parameters gamma_prob

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
%	      r =  crpi*(1-crr)*pinf +cry*(1-crr)*(y-yf)+crdy*(y-yf-y(-1)+yf(-1))+crr*r(-1)+ms  ;

% new variable created to re-introduce the monetary policy shock in optimal policy
% and avoid stochastic singularity under the estimation of the optimal policy model
		  dr=r-r(-1);
		  pinf_target=pinf+ms;
		  
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

observables dy dc dinve dw pinfobs robs labobs

planner_objective{discount = 0.99,commitment=gamma_prob} -.5*(1*pinf_target^2+wy*y^2+wr*dr^2);	

steady_state_model
	dy=ctrend;
	dc=ctrend;
	dinve=ctrend;
	dw=ctrend;
	pinfobs=constepinf;
	robs = conster;
	labobs = constelab;

parameterization
	% fixed parameters
	crhoas,1; 	 % this parameter does not enter the model
	crhols,    0.9928;  % this parameter does not enter the model
	ctou,.025;
	clandaw,1.5;
	cg,0.18;
	curvp,10;
	curvw,10;
	
	% estimated parameters initialisation
	sig_a,   0.4618,0.01,2,inv_gamma_pdf(0.9);
	sig_b,   0.18185,0.01,2,inv_gamma_pdf(0.9);
	sig_g,   0.6090,0.01,2,inv_gamma_pdf(0.9);
	sig_qs,  0.46017,0.01,2,inv_gamma_pdf(0.9);
	sig_m,   0.2397,0.01,2,inv_gamma_pdf(0.9);
	sig_pinf,0.1455,0.01,2,inv_gamma_pdf(0.9);
	sig_w,   0.2089,0.01,2,inv_gamma_pdf(0.9);

	calfa,   .24, .2,.4,normal_pdf(.9);
	csigma,   1.5,.75,2.25,normal_pdf(0.9);
	cfc,      1.5, 1, 1.5,normal_pdf(.9);
	cgy,      0.51,.1,1.5,normal_pdf(.9);
	
	csadjcost, 6.0144,1,7,normal_pdf(0.9);
	chabb,     0.6361,0.3,0.7,beta_pdf(.9);    
	cprobw,    0.8087,0.3,0.7,beta_pdf(.9);
	csigl,     1.9423,.5,3.5,normal_pdf(.9);
	cprobp,    0.6,0.3,0.7,beta_pdf(.9);
	cindw,     0.3243,0.3,0.7,beta_pdf(.9);
	cindp,     0.47,0.3,0.7,beta_pdf(.9);
	czcap,     0.2696,0.3,0.7,beta_pdf(.9);
	
%	% Taylor-rule monetary policy
%	crpi,      1.488,1,2,normal_pdf(.9);
%	crr,       0.8762,0.3,0.7,beta_pdf(.9);
%	cry,       0.0593,0.025,0.225,normal_pdf(.9);
%	crdy,      0.2347,0.025,0.225,normal_pdf(.9);

	% optimal policy
	wy, .6,.425,.825,gamma_pdf(.9);
	wr, .1,.2,.4,gamma_pdf(.9);
	
	crhoa,     0.9977,0.2,0.8,beta_pdf(.9);
	crhob,     0.5799,0.3,0.7,beta_pdf(.9);
	crhog,     0.9957,0.2,0.8,beta_pdf(.9);
	crhoqs,    0.7165,0.3,0.7,beta_pdf(.9);
	crhoms,    .3,0.3,0.7,beta_pdf(.9);
	crhopinf,  0.8,0.3,0.7,beta_pdf(.9);
	crhow,     0.8,0.3,0.7,beta_pdf(.9);
	cmap ,     0.7,0.3,0.7,beta_pdf(.9);
	cmaw  ,    0.7,0.3,0.7,beta_pdf(.9);
	
	% derived from steady state
	constebeta, 0.7420,.05,.8,gamma_pdf(.9);
	
	ctrend,     0.3982,.2,.6,normal_pdf(.9);
	constepinf, .7,.425,.825,gamma_pdf(.9);
	constelab,  1.2918,-4,4,normal_pdf(.9);

	gamma_prob,.5,0.3,0.7,beta_pdf(.9);

