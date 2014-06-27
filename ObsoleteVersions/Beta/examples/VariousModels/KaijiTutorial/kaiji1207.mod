var C, K, H, R, A, B, D, PSI;

varexo EPS_A, EPS_B, EPS_D, EPS_PSI;

parameters alpha, bss, dss, psiss, rho_psi, rho_a, rho_b, rho_d, sig_d, sig_a, sig_b, sig_psi;

varobs C, H, K;

model
	1/C = B*(1+R{+1}-D{+1})/C{+1};

	C + K = A*K{-1}^alpha*H^(1-alpha)+(1-D)*K{-1};

	(1-alpha)*A*(K{-1}/H)^alpha/C = PSI/(1-H);

	R = alpha*A*(K{-1}/H)^(alpha-1);

	log(A) = rho_a*log(A{-1})+sig_a*EPS_A;

	log(B/bss) = rho_b*log(B{-1}/bss)+sig_b*EPS_B;

	log(D/dss) = rho_d*log(D{-1}/dss)+sig_d*EPS_D;

	log(PSI/psiss) = rho_psi*log(PSI{-1}/psiss)+sig_psi*EPS_PSI;

end

steady_state_model
	A  = 1;
	B  = bss;
	D  = dss;
	PSI=psiss;
	C  = 0.6;
	K  = 2.0;
	H  = .31;
	R  = .04;
end

parameterization
	alpha,   .36;
	bss,     .96;
	dss,      .1;
	psiss,	1.45;
	rho_psi, .75, 0, 1;
	rho_a,   .96, 0, 1;
	rho_b,   .75, 0, 1;
	rho_d,   .75, 0.25, .75, beta_pdf, .9;
	sig_d,   .01, 0, 0.5;
	sig_a,   .01, 0, 0.5;
	sig_b,   .01, 0, 0.5;
	sig_psi, .01, 0, 0.5;
end
	