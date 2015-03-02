endogenous Y C PAI R G X Z A THETA

exogenous E_A, E_THETA, E_Z, E_R

parameters phi rho_a sig_a rho_theta sig_theta sig_z beta eta sig_r	rho_x rho_g rho_pai	rho_r
 		   pai_ss a_ss theta_ss	z_ss

model
	Y = C + phi/2*(PAI/steady_state(PAI)-1)^2*Y;

	log(A) = (1-rho_a)*log(a_ss)+rho_a*log(A{-1})+sig_a*E_A;

	A/C = beta*R*A{+1}/C{+1}*1/Z{+1}*1/PAI{+1};

	log(THETA) = (1-rho_theta)*log(theta_ss)+rho_theta*log(THETA{-1})+sig_theta*E_THETA;

	log(Z) = log(z_ss)+sig_z*E_Z;

	0 = 1-THETA+THETA*(C/A)*Y^(eta-1)-phi*(PAI/steady_state(PAI)-1)*PAI/steady_state(PAI)+
		beta*phi*A{+1}/A*C/C{+1}*(PAI{+1}/steady_state(PAI)-1)*PAI{+1}/steady_state(PAI)*Y{+1}/Y;

	G = Y/Y{-1}*Z;

	X = Y/A^(1/eta);

	R/steady_state(R) = (R{-1}/steady_state(R))^rho_r*
						(PAI/steady_state(PAI))^rho_pai*
						(G/steady_state(G))^rho_g*
						(X/steady_state(X))^rho_x*
						exp(sig_r*E_R);

observables G PAI R

steady_state_model
	Z = z_ss;
	A = a_ss;
	THETA = theta_ss;
	Y = (A*(THETA-1)/THETA)^(1/eta);
	C = Y;
	PAI	= pai_ss;
	R = PAI*Z/beta;
	G = Z;
	X = Y/A^(1/eta); % (A*(THETA-1)/THETA)^(1/eta);

parameterization
	% fixed parameters
	z_ss      ,1.0048;
	pai_ss    ,1.0086;
	beta      ,0.9900;
	phi       ,0.1000;
	rho_r     ,1.0000;
	theta_ss  ,5.0000;
	a_ss      ,1.0000;
	eta       ,1/0.0617;
	
	% omega     ,0.0617,0,1;
	% alpha_x   ,0.0836,0,1; 
	% alpha_pai ,0.0001,0,1; 
	rho_pai   ,0.3597,0,2; 
	rho_g     ,0.2536,0,2; 
	rho_x     ,0.0347,0,2; 
	rho_a     ,0.9470,0,2; 
	rho_theta ,0.9625,0,2; 
	sig_a     ,0.0405,0,2; 
	sig_theta ,0.0012,0,2; 
	sig_z     ,0.0109,0,2; 
	sig_r     ,0.0031,0,2; 
	