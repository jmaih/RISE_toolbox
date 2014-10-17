function flag=is_stable_system(obj,fast)
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

if nargin<2
    fast=true;
end
% this function checks that the solved system is
% stable in the sense that its covariance matrix is
% bounded.
% References:
% 1- Vijay Gupta, Richard Murray and Babak Hassibi (2003):
% "On the Control of Jump Linear Markov Systems with Markov State Estimation"
% 2- O.L.V. Costa, M.D. Fragoso and R.P. Marques (2004):
% "Discrete-Time Markov Jump Linear Systems"
if isempty(obj)
    flag=struct('stability_criterion',1.000001);
    return
end

nobj=numel(obj);
if nobj>1
    flag=nan(1,nobj);
    for iobj=1:nobj
        flag(iobj)=is_stable_system(obj(iobj));
    end
    return
end
% this is hard-coded for the moment. We don't want to check stability under
% estimation if the problem is too big...
ms_stability_check_threshold=5000;

% trim down in case some matrices are the same and keep only the states
[T,Q,n,h]=problem_reduction(obj.solution.Tz,obj.solution.transition_matrices.Q);

max_eig=0;
if ~isempty(T{1})
    % update n right here right now
    if h==1
        max_eig=max(abs(eig(full(T{1}))));
    else
        n2=n^2;
        if ~(obj.estimation_under_way && ...
                h*n2>ms_stability_check_threshold)
            % do not check stability under estimation if matrix is too big. I
            % should probably decrease the threshold coz what is expensive is
            % the calculation of the kronecker products. Moreover, I could put
            % the threshold as an option...
            C=gupta_murray_hassibi();
            if obj.options.debug
                C2=costa_fragoso_marques();
                test=max(max(abs(C-C2)));
                disp([mfilename,':: test of accuracy of computation of Mean Square Stability criterion ',num2str(test)])
            end
            max_eig=max(abs(eig(C)));
        end
    end
end

flag=max_eig<=obj.options.stability_criterion;%1+1e-12

% tic
% flag2=test_is_stable_system_old(obj.T,obj.Q);
% toc % verdict: at least 4 times faster than old
% if ~isequal(flag2,flag)
%     keyboard
% end

    function C=gupta_murray_hassibi()
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
    function C=costa_fragoso_marques()
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

    function [T_red,Q_red,n,h_red]=problem_reduction(T,Q)
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
        ov=obj.order_var.after_solve;
        pb=obj.locations.after_solve.t.pb;
        zpb=obj.locations.after_solve.z.pb;
        %================= trimming preparation=================%
        for i_state=1:h_red
            T_red{i_state}=T_red{i_state}(ov(pb),zpb);
            % target the variables that are known to be unstable a priori
            if ~obj.is_stationary_model
                if i_state==1
                    balanced_growth=abs(full(obj.solution.bgp{1}))>1e-4;
                    balanced_growth=balanced_growth(ov(pb));
                end
                T_red{i_state}=T_red{i_state}(~balanced_growth,~balanced_growth);
            end
        end
        %=======================================================%
        % update n
        n=size(T_red{1},1);
    end
end
