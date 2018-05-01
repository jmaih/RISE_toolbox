function [T_red,Q_red,n,h_red]=problem_reduction(obj)
% PROBLEM_REDUCTION - prepares the solution for checking Mean Square
% Stability of the system
%
% ::
%
%
%   [T,Q,n,h]=problem_reduction(obj)
%
% Args:
%
%    - **obj** [dsge|rise]: solved model object
%
% Returns:
%    :
%
%    - **T** [1 x h cell]: Companion form for different regimes. Each cell
%    contains an n x n matrix
%
%    - **Q** [h x h matrix]: Transition matrix
%
%    - **n** [integer]: number of endogenous variable
%
%    - **h** [integer]: number of regimes
%
% Note:
%
% Example:
%
%    See also: SVAR/PROBLEM_REDUCTION

T=obj.solution.Tz;
Q=obj.solution.transition_matrices.Q;
% do this when we have more than one markov chain
h=numel(T);
n=size(T{1},1);
if n==0
    warning([mfilename,':: the model has not been solved']) %#ok<WNTAG>
    T_red=T;
    Q_red=Q;
    h_red=h;
    return
end
regimes=1:h;
for istate=1:h
    if istate>1
        alone=true;
        for jstate=1:istate-1
            if max(max(abs(T{jstate}-T{istate})))<1e-10
                alone=false;
                regimes(istate)=regimes(jstate);
            end
            if ~alone
                break
            end
        end
    end
end

regimes_red=unique(regimes);
h_red=numel(regimes_red);
if isequal(regimes_red,regimes) % problem cannot be reduced
    Q_red=Q;
else
    Q_red=nan(h_red);
    for istate=1:h_red
        tmp=sum(Q(:,regimes==regimes_red(istate)),2);
        Q_red(:,istate)=tmp(regimes_red);
    end
end
T_red=T(regimes_red);
% Now focus on the states only
ov=obj.order_var;
pb=obj.locations.after_solve.t.pb;
zpb=obj.locations.after_solve.z.pb;
%================= trimming preparation=================%
for i_state=1:h_red
    T_red{i_state}=T_red{i_state}(ov(pb),zpb);
end
%=======================================================%
% update n
n=size(T_red{1},1);
end