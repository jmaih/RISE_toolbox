/*
 * This file replicates the estimation of the cash in advance model described
 * Frank Schorfheide (2000): "Loss function-based evaluation of DSGE models",
 * Journal of Applied Econometrics, 15(6), 645-670.
 *
 * The data are in file "fsdat_simul.m", and have been artificially generated.
 * They are therefore different from the original dataset used by Schorfheide.
 *
 * The equations are taken from J. Nason and T. Cogley (1994): "Testing the
 * implications of long-run neutrality for monetary business cycle models",
 * Journal of Applied Econometrics, 9, S37-S70.
 * Note that there is an initial minus sign missing in equation (A1), p. S63.
 *
 * This implementation was written by Michel Juillard. Please note that the
 * following copyright notice only applies to this Dynare implementation of the
 * model.
 */

/*
 * Copyright (C) 2004-2010 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

var m P c e W R k d n l gy_obs gp_obs y dA;
varexo e_a e_m;

parameters alp bet gam mst rho psi del;

//// starting point for estimation
//alp=  0.3560;
//bet=  0.9930;
//gam=  0.0085;
//mst=  1.0002;
//rho=  0.1290;
//psi=  0.6500;
//del=  0.0100;

// baseline calibration
alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;

model;
dA = exp(gam+e_a);
log(m) = (1-rho)*log(mst) + rho*log(m(-1))+e_m;
-P/(c(+1)*P(+1)*m)+bet*P(+1)*(alp*exp(-alp*(gam+log(e(+1))))*k^(alp-1)*n(+1)^(1-alp)+(1-del)*exp(-(gam+log(e(+1)))))/(c(+2)*P(+2)*m(+1))=0;
W = l/n;
-(psi/(1-psi))*(c*P/(1-n))+l/n = 0;
R = P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(-alp)/W;
1/(c*P)-bet*P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)/(m*l*c(+1)*P(+1)) = 0;
c+k = exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)+(1-del)*exp(-(gam+e_a))*k(-1);
P*c = m;
m-1+d = l;
e = exp(e_a);
y = k(-1)^alp*n^(1-alp)*exp(-alp*(gam+e_a));
gy_obs = dA*y/y(-1);
gp_obs = (P/P(-1))*m(-1)/dA;
end;

initval;
	k = 6;
	m = mst;
	P = 2.25;
	c = 0.45;
	e = 1;
	W = 4;
	R = 1.02;
	d = 0.85;
	n = 0.19;
	l = 0.86;
	y = 0.6;
	gy_obs = exp(gam);
	gp_obs = exp(-gam);
	dA = exp(gam);
end;

shocks;
//	// parameters used for estimation
//	var e_a; stderr 0.035449;
//	var e_m; stderr 0.008862;
	// baseline calibration
	var e_a; stderr  0.014;//
	var e_m; stderr  0.005;//
end;

steady;

check;

estimated_params;
	alp       , beta_pdf     , 0.356000, 0.02;
	bet       , beta_pdf     , 0.993000, 0.002;
	gam       , normal_pdf   , 0.008500, 0.003;
	mst       , normal_pdf   , 1.000200, 0.007;
	rho       , beta_pdf     , 0.129000, 0.223;
	psi       , beta_pdf     , 0.650000, 0.05;
	del       , beta_pdf     , 0.010000, 0.005;
	stderr e_a, inv_gamma_pdf, 0.035449, inf;
	stderr e_m, inv_gamma_pdf, 0.008862, inf;
end;

varobs gp_obs gy_obs;

stoch_simul(order=1)m P c e W R k d n l ;

return

estimation(order=1, datafile=fsdat_simul, nobs=192, mh_replic=2000, mh_nblocks=2, mh_jscale=0.8, mode_compute=1);  //, loglinear


/*
 * The following lines were used to generate the data file. If you want to
 * generate another random data file, comment the "estimation" line and uncomment
 * the following lines.
 */

//stoch_simul(periods=200, order=1);
//datatomfile('fsdat_simul', char('gy_obs', 'gp_obs'));
