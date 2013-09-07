exogenous E_WBAR E_Y

endogenous Y RHO I PAI W P V WBAR DI

parameters lamb nu lamb_pai


model(linear)

	Y = 1.45*Y(-1)-0.47*Y(-2)-0.34*RHO(-1)+.01*E_Y;
	
	RHO = 40/41*RHO(+1)+1/41*(I-PAI(+1));
	
	P = 0.42*W + 0.31*W(-1)+ 0.19*W(-2)+ 0.08*W(-3);
	
	PAI = 4*(P-P(-1));
	
	V = 0.42*WBAR + 0.31*WBAR(-1) + 0.19*WBAR(-2) + 0.08*WBAR(-3);
	
	WBAR = W-P;
	
	DI = I-I(-1);
		
	WBAR = 0.42*V + 0.31*V(+1) + 0.19*V(+2) + 0.08*V(+3)
	       +0.002*(0.42*Y + 0.31*Y(+1) + 0.19*Y(+2) + 0.08*Y(+3)) + .01*E_WBAR;
		   
parameterization
	lamb,1    ;
	nu,.5     ;
	lamb_pai,1;

planner_objective{discount = 0.99,commitment=1} -.5*(lamb_pai*PAI^2+lamb*Y^2+nu*DI^2+0*I^2);

