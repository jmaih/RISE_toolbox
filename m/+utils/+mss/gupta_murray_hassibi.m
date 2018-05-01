function [flag,max_eig,C]=gupta_murray_hassibi(T,Q,crit,fast)
% GUPTA_MURRAY_HASSIBI -- checks for Mean Square Stability using the Gupta,
% Murray and Hassibi algorithm
%
% ::
%
%
%   C=GUPTA_MURRAY_HASSIBI(T,Q)
%
%   C=GUPTA_MURRAY_HASSIBI(T,Q,crit)
%
%   C=GUPTA_MURRAY_HASSIBI(T,Q,crit,fast)
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
%    - **fast** [false|{true}]: use the fastest algorithm
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
%    See also: COSTA_FRAGOSO_MARQUES

% References:
% Vijay Gupta, Richard Murray and Babak Hassibi (2003):
% "On the Control of Jump Linear Markov Systems with Markov State Estimation"

if nargin<4
    fast=[];
    if nargin<3
        crit=[];
    end
end
if isempty(crit)
    crit=1+1e-12;
end
if isempty(fast)
    fast=true;
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
        n2=n^2;
        C=zeros(h*n2);
        TT=cell(1,h);
        for jstate=1:h
            jspan=(jstate-1)*n2+1:jstate*n2;
            sameas=[];
            for kstate=1:jstate-1
                if max(max(abs(T{jstate}-T{kstate})))<1e-9
                    sameas=kstate;
                    break
                end
            end
            if isempty(sameas)
                TT{jstate}=kron(T{jstate},T{jstate});
            else
                TT{jstate}=TT{sameas};
            end
            if fast
                C(:,jspan)=kron(Q(:,jstate),TT{jstate});
            else
                for istate=1:h
                    ispan=(istate-1)*n2+1:istate*n2;
                    C(ispan,jspan)=Q(istate,jstate)*TT{jstate};
                end
            end
        end
    end
end