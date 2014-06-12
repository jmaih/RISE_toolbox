function Phi=draw_Phi(Lambdat,priors)
%%  Code to draw (Q|B_t)
% Draw Phi, the covariance of Lambdat (from iWishart)
% priors: priors for mean and variance of Phi
% B: time varying B
%%
Phi_prmean=priors.Phimean;
Phi_prvar=priors.Phivar;
[t,M]=size(Lambdat);
Sigt=ln(Lambdat);
% Get first differences of Sigtdraw to compute the SSE
Sigttemp = Sigt(2:t,:) - Sigt(1:t-1,:);
sse_2 = zeros(M,M);
for i = 1:t-1
    sse_2 = sse_2 + Sigttemp(i,:)'*Sigttemp(i,:);
end
Phiinv = inv(sse_2 + Phi_prmean);
Phiinvdraw = wish(Phiinv,t+Phi_prvar);
Phi = inv(Phiinvdraw);  % this is a draw from W