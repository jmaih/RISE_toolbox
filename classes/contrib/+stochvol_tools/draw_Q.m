function Q=draw_Q(B,priors)
%%  Code to draw (Q|B_t)
% Draw Q, the covariance of Bt (from iWishart)
% priors: priors for mean and variance of Q
% B: time varying B
%%
Q_prmean=priors.Qmean;
Q_prvar=priors.Qvar;
K=size(Q_prmean,1);
% Take the SSE in the state equation of B(t)
Btemp = B(2:t,:) - B(1:t-1,:);
sse_2 = zeros(K,K);
for i = 1:t-1
    sse_2 = sse_2 + Btemp(i,:)'*Btemp(i,:);
end

% ...and subsequently draw Q, the covariance matrix of B(t)
Qinv = inv(sse_2 + Q_prmean);
Qinvdraw = wish(Qinv,t+Q_prvar);
Q = inv(Qinvdraw);  % this is a draw from Q