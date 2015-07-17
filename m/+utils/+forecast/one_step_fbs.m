function [y1,is_active_shock,retcode,shocks]=one_step_fbs(T,y0,ss,xloc,sig,...
    shocks,order,compl,cond_shocks_id)

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

narginchk(7,9)

n=nargin;
[nx,kplus1]=size(shocks);

myzero=-sqrt(eps);


y1k=y0(:,ones(1,kplus1+1));
if n==7
    compl=[];
    nrows=0;
else
    nrows=size(compl(y0.y),1);
end

% compute the first forecast and detect any violation: quick exit if none
%------------------------------------------------------------------------
check_first_only=true;
[~,first_viol]=multi_step_shooting(shocks,true);

y1=y1k(2);
retcode=0;
is_active_shock=any(abs(shocks)>sqrt(eps),1);
if ~isempty(first_viol)
    debug=false;
    if debug
        fsolve_options=struct('Display','iter','TolFun',1e-12,'TolX',1e-12);
    else
        fsolve_options=struct('Display','none','TolFun',1e-12,'TolX',1e-12);
    end
    % now check the entire stretch
    %------------------------------
    check_first_only=false;
    is_active_shock(first_viol)=true;
    if islogical(cond_shocks_id)
        cond_shocks_id=find(cond_shocks_id);
    end
    n_cond_shocks=numel(cond_shocks_id);
    if n_cond_shocks~=nrows
        error('Singularity in the problem of matching restrictions')
    end
    done=false;
    while ~done
        % find the shocks that make the violations go away and the
        % constraints bind exactly
        %------------------------------------------------------------------
        shocks0=shocks;
        shocks0(cond_shocks_id,is_active_shock)=nan;
        e_id=isnan(shocks0);
        ee0=zeros(n_cond_shocks*sum(is_active_shock),1);
        if debug
            disp(is_active_shock)
        end
        [ee,fval,exitflag]=fsolve(@multi_complementarity,ee0,fsolve_options);
        exitflag=utils.optim.exitflag(exitflag,ee,max(abs(fval(:))));
        if exitflag~=1
            retcode=701;
            return
        end
        shocks0(e_id)=ee;
        % given those shocks, check that all future steps do not violate
        % the constraints
        %------------------------------------------------------------------
        [~,first_viol]=multi_step_shooting(shocks0,true);
        
        test_passed=isempty(first_viol);
        % if the test is passed we are done
        %----------------------------------
        if test_passed
            done=true;
            if all(is_active_shock)
                % the last guys is also constrained. Do one more step to see if
                % there is a lift up. If no lift up we are not out of the bush
                %--------------------------------------------------------------
                y1_last_uncond=utils.forecast.one_step_engine(T,y1k(end),ss,...
                    xloc,sig,zeros(size(shocks0)),order);
                if ~all(compl(y1_last_uncond.y)>myzero)
                    retcode=701;
                    return
                end
            end
        else
            is_active_shock(first_viol)=true;
        end
    end
    y1=y1k(2);
    shocks=shocks0;
end

    function viol=multi_complementarity(ee)
        shocks0(e_id)=ee;
        [viol,first_viol]=multi_step_shooting(shocks0);
        viol=viol(:,is_active_shock);
        viol=viol(:);
    end

    function [viol,first_viol]=multi_step_shooting(shocks,stop_at_viol)
        if nargin<2
            stop_at_viol=false;
        end
        viol=zeros(nrows,kplus1);
        first_viol=[];
        for icol_=2:kplus1+1 %<--size(y,2)
            shocks_i=[shocks(:,icol_-1:end),zeros(nx,icol_-2)];
            y1k(icol_)=utils.forecast.one_step_engine(T,y1k(icol_-1),ss,...
                xloc,sig,shocks_i,order);
            if isempty(compl)
                % if no restrictions, exit after the first forecast to
                % save time 
                return
            end
            viol(:,icol_-1)=compl(y1k(icol_).y);
            if isempty(first_viol) && ~all(compl(y1k(icol_).y)>myzero)
                first_viol=icol_-1;
                if stop_at_viol
                    break
                end
            end
            if icol_==2 && check_first_only
                return
            end
        end
    end
end