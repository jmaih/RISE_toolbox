function [flag,max_eig,C]=costa_fragoso_marques(T,Q,crit)
% COSTA_FRAGOSO_MARQUES -- checks for Mean Square Stability using the
% Costa, Fragoso and Marques algorithm.
%
% ::
%
%
%   C=COSTA_FRAGOSO_MARQUES(T,Q)
%
%   C=COSTA_FRAGOSO_MARQUES(T,Q,crit)
%
% Args:
%
%    - **T** [1 x h cell array]: Autoregressive part of the solution in each
%    regime
%
%    - **Q** [h x h matrix]: Transition matrix
%
%    - **crit** [numeric]: tolerance criterion
%
% Returns:
%    :
%
%    - **flag** [true|false]: result of the investigation on whether there is
%    MSS or not.
%
%    - **max_eig** [numeric]: maximum eigenvalue
%
%    - **C** [matrix]: matrix in the criterion for checking stability
%
% Note:
%
% Example:
%
%    See also: GUPTA_MURRAY_HASSIBI

% References:
% O.L.V. Costa, M.D. Fragoso and R.P. Marques (2004):
% "Discrete-Time Markov Jump Linear Systems"
if nargin<3
    crit=[];
end
if isempty(crit)
    crit=1+1e-12;
end

n=size(T{1},1);
h=size(Q,1);

if h==1
    C=do_one_regime();
else
    C=do_multiple_regimes();
end

max_eig=max(abs(eig(full(C))));

flag=max_eig<=crit;

    function C=do_one_regime()
        C=T{1};
    end

    function C=do_multiple_regimes()
        n=size(T{1},1);
        h=size(Q,1);
        n2=n^2;
        D=zeros(h*n2);
        for istate=1:h
            ispan=(istate-1)*n2+1:istate*n2;
            sameas=[];
            if istate>1
                for ii=1:istate-1
                    if max(max(abs(T{istate}-T{ii})))<1e-9
                        sameas=ii;
                        break
                    end
                end
            end
            notdone=istate==1||isempty(sameas);
            if notdone
                tmp=kron(T{istate},T{istate});
            else
                iter_same=(sameas-1)*n2+1:sameas*n2;
                tmp=D(iter_same,iter_same);
            end
            D(ispan,ispan)=tmp;
        end
        C=kron(Q,eye(n2))*D;
    end
end
