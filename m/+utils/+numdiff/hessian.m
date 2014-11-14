function H = hessian(Objective,xparam,varargin)
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

tol      = eps.^(1/4);
diagonly = false;

Objective=fcnchk(Objective,length(varargin));

fx = Objective(xparam,varargin{:});

npar = size(xparam,1);

% Compute the stepsize (h)
h = tol.*max(abs(xparam),1);
xh = xparam+h;
h = xh-xparam;
ee = sparse(1:npar,1:npar,h,npar,npar);

hh=h*h';

% Compute forward and backward steps
f1 = zeros(npar,1);
f0 = zeros(npar,1);
Hdiag = zeros(npar,1);
theLoopBody1=@loop_body_diagonal;
theLoopBody2=@loop_body_cross;
some_workers=false;
%try
%	some_workers=~isempty(gcp('nocreate'));
%catch
%	some_workers=matlabpool('size');
%end
if license('checkout','Distrib_Computing_Toolbox') && some_workers
    parfor ii=1:npar
        iter=ii;
        [f1(ii),f0(ii),Hdiag(ii)]=theLoopBody1(iter);
    end
    H=diag(Hdiag);
    if ~diagonly
        parfor ii=1:npar
            H(ii,:)=theLoopBody2(H(ii,:),ii);
        end
    end
else
    for ii=1:npar
        [f1(ii),f0(ii),Hdiag(ii)]=theLoopBody1(ii);
    end
    H=diag(Hdiag);
    if ~diagonly
        for ii=1:npar
            H(ii,:)=theLoopBody2(H(ii,:),ii);
        end
    end
end
if ~diagonly
    H=H+triu(H,1)';
    H=.5*(H+H');
end

    function [f1i,f0i,Hii]=loop_body_diagonal(ii)
        x1 = xparam+ee(:,ii);
        f1i = Objective(x1,varargin{:});
        x0 = xparam-ee(:,ii);
        f0i = Objective(x0,varargin{:});
        Hii = (f1i+f0i-2*fx)./hh(ii,ii);
    end
    function H=loop_body_cross(H,ii)
        for jj=ii+1:npar
            xcross =  xparam+ee(:,ii)-ee(:,jj);
            fxx=Objective(xcross,varargin{:});
            H(jj) = (f1(ii)+f0(jj)-fx-fxx)./hh(ii,jj);
        end
    end
end