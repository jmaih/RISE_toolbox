function varargout=dirichlet(varargin)

hyperparameter_mode=nargin>0;
if hyperparameter_mode
    [m,s]=hyperparameters_2_moments(varargin{:});
    a=varargin{1};
    moments=struct('mean',m,'sd',s);
    % Let b be the complementary to the marginal
    b=sum(a)-a;
    fval=0;
    [~,space]=hyperparameters(a);
    varargout={a,b,moments,fval,space};
else
    lpdfn=@logdens_dirichlet;
    cdfn=@(varargin)[];
    icdfn=@(varargin)[];
    rndfn=@draw_dirichlet;
    m2h=@moments_2_hyperparameters;
    h2m=@hyperparameters_2_moments;
    varargout={lpdfn,cdfn,icdfn,rndfn,m2h,h2m};
end

end

function lpdf=logdens_dirichlet(theta,a,~,~,d)
if nargin<5
    d=1;
end
% account for truncation
%------------------------
theta=theta/d;
if all(theta>=0)
    lpdf=gammaln(sum(a))-sum(gammaln(a(:)))+sum((a(:)-1).*log(theta(:)));
else
    lpdf=-inf;
end
end

function draw=draw_dirichlet(a,~,n,~,d)
if nargin<5
    d=1;
    if nargin<3
        n=1;
    end
end
k=numel(a);
u=rand(k,n);
aa=a(:);
draw=gaminv(u,aa(:,ones(1,n)),1);
draw=bsxfun(@rdivide,draw,sum(draw,1));
% make small as necessary to account for a possible truncation
%---------------------------------------------------------------
draw=draw*d;
end

function [m,s]=hyperparameters_2_moments(varargin)
a=varargin{1};
if nargin==4
    d=varargin{4};
else
    d=1;
end
m=a(:)/sum(a);
s=sqrt(a(:).*(sum(a)-a(:))/(sum(a)^2*(sum(a)+1)));
% account for truncation
%-----------------------
m=d*m;
s=d*s;
end

function [a,b]=moments_2_hyperparameters(m,s,~,d)
if nargin<4
    d=1;
end
% account for truncations
%------------------------
m=m(:)/d;
s=s(:)/d;
b=nan(size(m));
% 1- find sum of hyperparameters
%--------------------------------
a0=m.*(1-m)./s.^2-1;
if any(abs(a0-a0(1))>1e-10)||a0(1)<=0
    error('wrong moments')
end
% 2- find the hyperparameters
%------------------------------
a=a0.*m;
end

function [violation,space]=hyperparameters(hyper)
% -inf<a<inf, b>0
space=[1e-12,inf];
if nargin==0||isempty(hyper)
    violation=space;
else
    violation=any(hyper(:)<space(:,1))||any(hyper(:)>space(:,2));
end
end