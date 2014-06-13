function H = hessian(Objective,xparam,varargin)
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
if license('checkout','Distrib_Computing_Toolbox') && matlabpool('size')
    parfor ii=1:npar
        iter=ii;
        [f1(ii),f0(ii),Hdiag(ii)]=theLoopBody1(iter);
    end
    H=diag(Hdiag);
    if ~diagonly
        parfor ii=1:npar
            H(ii,:)=theLoopBody2(H(ii,:),f1,f0,ii);
        end
    end
else
    for ii=1:npar
        [f1(ii),f0(ii),Hdiag(ii)]=theLoopBody1(ii);
    end
    H=diag(Hdiag);
    if ~diagonly
        for ii=1:npar
            H(ii,:)=theLoopBody2(H(ii,:),f1,f0,ii);
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
    function H=loop_body_cross(H,f1,f0,ii)
        for jj=ii+1:npar
            xcross =  xparam+ee(:,ii)-ee(:,jj);
            fxx=Objective(xcross,varargin{:});
            H(jj) = (f1(ii)+f0(ii)-fx-fxx)./hh(ii,jj);
        end
    end
end


% if exist('matlabpool.m','file') && matlabpool('size')>0
% %     disp([mfilename,':: using parallel code'])
%     diag_terms=nan(npar,1);
%     diag_locs=(0:npar-1)*npar+(1:npar);
%     hh_diag=hh(diag_locs);
%     parfor i=1:npar
%         x1 = xparam+ee(:,i);
%         f1(i) = Objective(x1,varargin{:}); %#ok<PFBNS>
%         x0 = xparam-ee(:,i);
%         f0(i) = Objective(x0,varargin{:});
%         diag_terms(i)= (f1(i)+f0(i)-2*fx)./hh_diag(i);
%     end
%     H(diag_locs)=diag_terms;
%    % Compute double steps
%     if ~diagonly
%         parfor i=1:npar
%             Hij=zeros(npar,1);
%             Hij(i)=diag_terms(i);
%             hh_i=hh(i,:);
%             for j=i+1:npar
%                 xcross =  xparam+ee(:,i)-ee(:,j); %#ok<PFBNS>
%                 fxx=Objective(xcross,varargin{:}); %#ok<PFBNS>
%                 Hij(j)=(f1(i)+f0(j)-fx-fxx)./hh_i(j); %#ok<PFBNS>
%             end
%             H(i,:)=Hij;
%         end
%         H=H+triu(H,1)';
%     end
% else
%     disp([mfilename,':: using serial code'])
