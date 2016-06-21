endogenous y r pai g x a theta c z

exogenous e_a, e_theta, e_z, e_r

parameters rho_a, a_ss sig_a rho_theta theta_ss sig_theta
z_ss sig_z pai_ss beta sig_r omega rho_pai rho_g rho_x psi rho_r
 

model
	# eta=1/omega;
	# phi=eta*(theta_ss-1)/psi;
	
	log(a)=(1-rho_a)*log(a_ss)+rho_a*log(a{-1})+sig_a*e_a;
	
	log(theta)=(1-rho_theta)*log(theta_ss)+rho_theta*log(theta{-1})+sig_theta*e_theta;
	
	log(z)=log(z_ss)+sig_z*e_z;
	
	y=c+phi/2*(pai/pai_ss-1)^2*y;
	
	a/c=beta*r*a{+1}/c{+1}*(1/z{+1})*(1/pai{+1});
	
	0=1-theta+theta*(c/a)*y^(eta-1)-phi*(pai/pai_ss-1)*(pai/pai_ss)+beta*phi*a{+1}/a*c/c{+1}*(pai{+1}/pai_ss-1)*pai{+1}/pai_ss*(y{+1}/y);
	
	g=y/y{-1}*z;
	
	x=y/a^(1/eta);

	log(r/steady_state(r))=rho_r*log(r{-1}/$(r))+rho_pai*log(pai/pai_ss)+
		rho_g*log(g/steady_state(g))+rho_x*log(x/$(x))+sig_r*e_r;

	% steady_state(v) is the same as $(v)

parameterization
	% fixed parameters
	a_ss      ,1.0000;
	z_ss      ,1.0048;
	pai_ss    ,1.0086;
	beta      ,0.9900;
	psi       ,1.0000;
	rho_r     ,1.0000;
	theta_ss  ,6.0000;
	
	omega     ,0.0617;
%	alpha_x   ,0.0836; 
%	alpha_pai ,0.0000; 
	rho_pai   ,0.3597; 
	rho_g     ,0.2536; 
	rho_x     ,0.0347; 
	rho_a     ,0.9470; 
	rho_theta ,0.9625; 
	sig_a     ,0.0405; 
	sig_theta ,0.0012; 
	sig_z     ,0.0109; 
	sig_r     ,0.0031; 
