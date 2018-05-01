function H = hessian(Objective,xparam,options,varargin)
% hessian - computes the hessian of a scalar function
%
% ::
%
%
%   H = hessian(Objective,xparam)
%   H = hessian(Objective,xparam,varargin)
%
% Args:
%
%    - **Objective** [char|function handle]: function to differentiate
%
%    - **xparam** [vector]: Point at which the derivatives are taken
%
%    - **options** [struct]: structure containing the options for the
%      computation. These are:
%      - **tol** [numeric|{eps.^(1/4)}] : tolerance for the computation of the
%          steps
%      - **diagonly** [true|{false}] : if true, only the elements on the
%          diagonal are computed
%
%    - **varargin** []: further input arguments of **Objective**
%
% Returns:
%    :
%
%    - **H** [matrix]: hessian matrix
%
% Note:
%
%    - The tolerance level is hard-coded to be eps.^(1/4). This should
%      probably be an input
%
% Example:
%
%    See also:

if nargin<3
    options=[];
end
if isempty(options)
    options=struct();
end
if ~isfield(options,'tol')
    options.tol=[];
end
if ~isfield(options,'diagonly')
    options.diagonly=[];
end
defaults={
    'tol',eps.^(1/4),@(x)isnumeric(x)
    'diagonly',false,@(x)islogical(x)
    };

[tol,diagonly]=parse_arguments(defaults,...
    'tol',options.tol,'diagonly',options.diagonly);

if ischar(Objective)
    Objective=str2func(Objective);
end

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