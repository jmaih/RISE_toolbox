function PHI=vma(THETA,order,debug)
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

% Computes the Vector Moving Average representation of a reduced-form VAR

if nargin<3
    debug=false;
end

cell_mode=iscell(THETA);
if cell_mode
    nvar=size(THETA{1},1);
    nlags=numel(THETA);
    PHI=num2cell(zeros(1,order));
else
    [nvar,~,nlags]=size(THETA);
    PHI=zeros(nvar,nvar,order);
end

if debug
    fi_1=eye(nvar);
end

for io=1:order
    for ilag=1:min(nlags,io)
        jj=io-ilag;
        phi_j=1;
        if jj
            if cell_mode
                phi_j=PHI{jj};
            else
                phi_j=PHI(:,:,jj);
            end
        end
        if cell_mode
            PHI{io}=PHI{io}+THETA{ilag}*phi_j;
        else
            PHI(:,:,io)=PHI(:,:,io)+THETA(:,:,ilag)*phi_j;
        end
    end
    if debug
        if cell_mode
            fi_1=fi_1+PHI{io};
        else
            fi_1=fi_1+PHI(:,:,io);
        end
    end
end

if debug
    t_1=eye(nvar);
    for ilag=1:nlags
        t_1=t_1-THETA(:,:,ilag);
    end
    disp(max(max(abs(inv(fi_1)-t_1))))
end