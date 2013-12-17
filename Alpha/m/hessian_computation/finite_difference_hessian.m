function H = finite_difference_hessian(Objective,xparam,varargin)
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
H=zeros(npar,npar);

% Compute forward and backward steps
f1 = zeros(npar,1);
f0 = zeros(npar,1);
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
theLoopBody=@loop_body;
if matlabpool('size')>0
    parfor ii=1:npar
        theLoopBody()
    end
else
    for ii=1:npar
        theLoopBody()
    end
end
    if ~diagonly
        H=(H+H')./2;
    end
    function loop_body()
        x1 = xparam+ee(:,ii);
        f1(ii) = Objective(x1,varargin{:});
        x0 = xparam-ee(:,ii);
        f0(ii) = Objective(x0,varargin{:});
        H(ii,ii) = (f1(ii)+f0(ii)-2*fx)./hh(ii,ii);
        if ~diagonly
            % Compute double steps
            for jj=ii+1:npar
                if ii~=jj
                    xcross =  xparam+ee(:,ii)-ee(:,jj);
                    fxx=Objective(xcross,varargin{:});
                    H(ii,jj) = (f1(ii)+f0(jj)-fx-fxx)./hh(ii,jj);
                    H(jj,ii) =H(ii,jj);
                end
            end
        end
    end
end


