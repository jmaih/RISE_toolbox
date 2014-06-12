function Lambdat = draw_Lambda(y,B,A,Lambdatdraw,parameters,priors)
%% Code to draw (B_t|A,Lambda_t,Q)
% INPUT
% priors: priors for mean and variance of Sigma0=ln(Lambda0)
% data
% Zs = kron(ones(t,1),eye(M));
% B: time varying B
% Lambdatdraw: previous draw of time-varying volatilities
% Parameters
% Phi: variance/covariance matrix of shocks for ln(Lambdat)

persistent Zs smpl n
if isempty(Zs)
    [n,smpl]=size(y);
    Zs = kron(ones(smpl,1),eye(n));
end
%%
% call mixture of normals
[q_s,m_s,u2_s]=stochvol_tools.mixturevalues();

sigma_prmean=priors.sigma_0;
sigma_prvar=priors.sigma_var;
Phi=parameters.Phi;

% create yhat is the vector y(t) - Z x B(t) 
% create prediction error
%------------------------------------------
yhat = zeros(n,smpl);
for t = 1:smpl
    yhat(:,t) = y(:,t) - Z((t-1)*n+1:t*n,:)*B(t,:)';
end

% Then take squares
% of the resulting quantity (saved in matrix y2)
y2=(A\yhat).^2;

% Transform to the log-scale but also add the 'offset constant' to prevent
% the case where y2 is zero (in this case, its log would be -Infinity)
yss = log(y2 + 0.0000000001);

% Next draw statedraw (chi square approximation mixture component) conditional on Lambdadraw
% This is used to update at the next step the log-volatilities Sigtdraw
nmix=numel(m_s);
prw=nan(nmix,1);
statedraw=nan(smpl,n);
sig_u2_s=sqrt(u2_s);
for jj = 1:n
    for t = 1:smpl
        for kk = 1:nmix
             x=yss(t,jj);
             mu=Lambdatdraw(t,jj) + m_s(kk) - 1.2704;
            temp1=normpdf(x,mu,sig_u2_s(kk));%
%             temp1= (1/sqrt(2*pi*u2_s(kk)))*exp(-.5*(((yss(t,jj) - Lambdatdraw(t,jj) - m_s(kk) + 1.2704)^2)/u2_s(kk)));
            prw(kk) = q_s(kk)*temp1;
        end
        prw = prw./sum(prw);
        cprw = [0,cumsum(prw(:)')];
        cprw(end)=1;
        trand = rand;
        statedraw(t,jj)=find(cprw>trand,1,'first')-1;  % this is a draw of the mixture component index
    end
end

% In order to draw the log-volatilies, substract the mean and variance
% of the 7-component mixture of Normal approximation to the measurement
% error covariance
vart = zeros(smpl*n,n);
yss1 = zeros(smpl,n);
for t = 1:smpl
    for kk = 1:n
        imix = statedraw(t,kk);
        vart((t-1)*n+kk,kk) = u2_s(imix);
        yss1(t,kk) = yss(t,kk) - m_s(imix) + 1.2704;
    end
end

% Sigtdraw is a draw of the diagonal log-volatilies, which will give us SIGMA(t)
Sigtdraw = stochvol_tools.carter_kohn(yss1',Zs,vart,Phi,n,n,smpl,sigma_prmean,sigma_prvar);

% Draws in Sigtdraw are in logarithmic scale (log-volatilies). Create
% original standard deviations of the VAR covariance matrix
Lambdattemp = eye(n);
Lambdat = zeros(n*smpl,n);
for t = 1:smpl
    for kk = 1:n
        Lambdattemp(kk,kk) = exp(Sigtdraw(kk,t));
    end
    Lambdat((t-1)*n+1:t*n,:) = Lambdattemp;
end

    