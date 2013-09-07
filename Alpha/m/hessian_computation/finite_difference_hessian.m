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
    for i=1:npar
        x1 = xparam+ee(:,i);
        f1(i) = Objective(x1,varargin{:});
        x0 = xparam-ee(:,i);
        f0(i) = Objective(x0,varargin{:});
        H(i,i) = (f1(i)+f0(i)-2*fx)./hh(i,i);
    end
    % Compute double steps
    if ~diagonly
        for i=1:npar
            for j=i+1:npar
                if i~=j
                    xcross =  xparam+ee(:,i)-ee(:,j);
                    fxx=Objective(xcross,varargin{:});
                    H(i,j) = (f1(i)+f0(j)-fx-fxx)./hh(i,j);
                    H(j,i) =H(i,j);
                end
            end
        end
        H=(H+H')./2;
    end
% end


