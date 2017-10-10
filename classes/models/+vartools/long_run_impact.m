function Linf=long_run_impact(B,L0,nx)
% B: matrix of coefficients of the reduced-form VAR including the
% deterministic variables
%
% L0: Impact in period 0: square root of the covariance matrix of shocks
%
% nx: number of deterministic variables including the constant
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

n=size(B,1);

B=B(:,nx+1:end);

beta=0;

while ~isempty(B)
    
    beta_i=B(:,1:n);
    
    beta = beta + beta_i;
    
    B(:,1:n)=[];
    
end

Linf = (eye(n)-beta)\L0;

end