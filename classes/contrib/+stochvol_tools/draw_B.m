function Btdraw = draw_B(y,Z,A,Lambdat,Qdraw,priors)
%% Code to draw (B_t|A,Lambda_t,Q)
% INPUT
% priors: priors for mean and variance of B0
% data
% A covariance matrix
% Lambdat: time-varying volatilities
% Parameters
% Q: variance/covariance matrix of shocks for B

%%
B_0_prmean=priors.B_0;
B_0_prvar=priors.B_0_var;
K=size(B_0_prmean,2);
[t,M]=size(y);

lamb_cols=size(Lambdat,2);
Ht = zeros(M,M,lamb_cols);
for ii = 1:lamb_cols
    stem = Lambdat(:,ii);
    Hsd = A*diag(stem);
    Hdraw = Hsd*Hsd';
    Ht(:,:,ii) = Hdraw;  % H(t)
end
    
% Btdrawc is a draw of the mean VAR coefficients, B(t)
Btdrawc = stochvol_tools.carter_kohn(y',Z,Ht,Qdraw,K,M,t,B_0_prmean,B_0_prvar);

% Accept draw
Btdraw = Btdrawc';
    
% or use the code below to check for stationarity
% if stationary==1
%     %Now check for the polynomial roots to see if explosive
%     counter = false(1,t);
%     for ii = 1:t;
%         BBtempor = Btdrawc(:,ii);
%         BBtempor = reshape(BBtempor,M*p,M)';
%         ctemp1 = [BBtempor; eye(M*(p-1)) zeros(M*(p-1),M)];
%         counter(ii) =max(abs(eig(ctemp1)))>0.9999;
%     end
%     %if they have been rejected keep old draw, otherwise accept new draw
%     if sum(counter)==0
%         Btdraw = Btdrawc';
%     end
% end
