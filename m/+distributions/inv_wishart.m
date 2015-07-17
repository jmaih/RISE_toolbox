function varargout=inv_wishart(S,v,~,~)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
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


hyperparameter_mode=nargin>0;

if hyperparameter_mode
    if nargin<2
        error([mfilename,':: number of inputs should be either 0 or 2'])
    end
    k=size(S,1);
    [~,pp]=transpose(chol(S));
    if pp>0
        error([mfilename,':: matrix not positive definite'])
    end
    if v<=k+1
        error([mfilename,':: insufficient number of degrees of freedom'])
    end
    moments=struct('mean',1/(v-k-1)*S,'sd',[]);
    fval=0;
    space=[];
    varargout={S,v,moments,fval,space};
else
    lpdfn=@logdens_inverse_wishart;
    rndfn=@draw_inverse_wishart;
    cdfn=@(theta)[];
    icdfn=@(u)[];
    m2h=[];
    h2m=[];
    varargout={lpdfn,cdfn,icdfn,rndfn,m2h,h2m};
end

end

function lpdf=logdens_inverse_wishart(theta,S,v)
k=size(S,1);W=reshape(theta,k,k);
[~,pp]=chol(W);
if pp
	lpdf=-inf;
else
	G=prod(exp(gammaln(.5*(v+1-(1:k)))));
	lpdf=-log(2^(v*k)*pi^(0.25*k*(k-1))*G)...
	+0.5*v*log(det(S))...
	-0.5*(v+k+1)*log(det(W))...
	-0.5*trace(S/W);
end
end

function d=draw_inverse_wishart(S,v)
k=size(S,1);
if v<k
    error([mfilename,':: the degrees of freedom should exceed the dimension of S'])
end
mu=0; % <--- zeros(k,v)
C=chol(S,'lower');
X=mu+C*randn(k,v);
d=inv(X*X');
end

