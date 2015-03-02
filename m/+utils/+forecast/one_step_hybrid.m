function [y1,is_active_shock,retcode]=one_step_hybrid(T,y0,ss,xloc,sig,...
    shocks,order,compl,cond_shocks_id,shoot)

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

narginchk(7,10)

n=nargin;
if n<10
    shoot=[];
    if n<9
        cond_shocks_id=[];
        if nargin<8
            compl=[];
        end
    end
end
if isempty(compl)
    compl=@(x)1;
end
if isempty(shoot)
    shoot=false;
end
[nx,kplus1]=size(shocks);
if isempty(cond_shocks_id)
    cond_shocks_id=1:nx;
end

myzero=-sqrt(eps);
badvector=@(x)any(compl(x)<myzero);
debug=true;
if debug
    fsolve_options=struct('Display','iter','TolFun',1e-12,'TolX',1e-12);
else
    fsolve_options=struct('Display','none','TolFun',1e-12,'TolX',1e-12);
end

y1k=y0(:,ones(1,kplus1+1));
nrows=size(compl(y0.y),1);
% compute all the forecasts and detect the first violation
%---------------------------------------------------------
if shoot
    up_to=inf;
else
    % at the beginning check only the first guy
    %------------------------------------------
    up_to=1;
end
[~,first_viol]=multi_step_shooting(shocks,true);

% try a quick exit
%------------------
y1=y1k(2);
retcode=0;
is_active_shock=false(1,kplus1);
if ~isempty(first_viol)
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
                if badvector(y1_last_uncond.y)
                    retcode=701;
                    return
                end
            end
        else
            is_active_shock(first_viol)=true;
            % update the number of guys to check
            if ~shoot
                up_to=up_to+1;
            end
        end
    end
    y1=y1k(2);
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
        for icol_=1:kplus1 %<--size(y,2)
            shocks_i=[shocks(:,icol_:end),zeros(nx,icol_-1)];
            y1k(icol_+1)=utils.forecast.one_step_engine(T,y1k(icol_),ss,...
                xloc,sig,shocks_i,order);
            viol(:,icol_)=compl(y1k(icol_+1).y);
            if isempty(first_viol) && badvector(y1k(icol_+1).y)
                first_viol=icol_;
                if stop_at_viol
                    break
                end
            elseif icol_==up_to
                break
            end
        end
    end

end