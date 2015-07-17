% Woodford(2003) model as presented in Seonghoon Cho and Antonio Moreno (2011): 
% "The forward method as a solution refinement in rational expectations Models"
% Journal of Economics Dynamics and Control 35 (2011) 257-272
%
% The rise flag "multiple" when set to true gives the model with multiple equilibria

endogenous PAI, Y, R, I

exogenous EPS

parameters delta kappa mu, beta, phi, rho, lambda

model
	PAI = delta*PAI{+1}+kappa*Y;

	Y = mu*Y{+1}+(1-mu)*Y{-1}-phi*(I-PAI{+1}-R);

	I = beta*PAI{+1}+lambda*Y;

	R = rho*R{-1} + EPS;

parameterization
	delta, 0.99;
	kappa,0.3;
	mu, .55;
	@#if multiple
		beta,0.95;
	@#else
		beta,1.5;
	@#end
	phi, 1;
	rho,.8;
	lambda,.1;