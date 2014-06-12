/*
 * This file is a modified version of 'fs2000.mod'.
 *
 * The difference is that, here, the equations are written in non-stationary form,
 * and Dynare automatically does the detrending.
 *
 * Also note that "m" and "dA" in 'fs2000.mod' are here called "gM" and "gA"
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

var gM gA;
trend_var(growth_factor=gA) A;
trend_var(growth_factor=gM) M;
var(deflator=A) k c y;
var(deflator=M(-1)/A) P;
var(deflator=M(-1)) W l d;
var R n;
varexo e_a e_m;

parameters alp bet gam mst rho psi del;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;

model;
gA = exp(gam+e_a);
log(gM) = (1-rho)*log(mst) + rho*log(gM(-1))+e_m;
c+k = k(-1)^alp*(A*n)^(1-alp)+(1-del)*k(-1);
P*c = M;
P/(c(+1)*P(+1))=bet*P(+1)*(alp*k^(alp-1)*(A(+1)*n(+1))^(1-alp)+(1-del))/(c(+2)*P(+2));
(psi/(1-psi))*(c*P/(1-n))=W;
R = P*(1-alp)*k(-1)^alp*A^(1-alp)*n^(-alp)/W;
W = l/n;
M-M(-1)+d = l;
1/(c*P)=bet*R/(c(+1)*P(+1));
y = k(-1)^alp*(A*n)^(1-alp);
end;

initval;
k = 6;
gM = mst;
P = 2.25;
c = 0.45;
W = 4;
R = 1.02;
d = 0.85;
n = 0.19;
l = 0.86;
y = 0.6;
gA = exp(gam);
end;

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
end;

steady;

check;

stoch_simul;
