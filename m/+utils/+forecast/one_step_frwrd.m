function [y1,is_active_shock,retcode]=one_step_frwrd(T,y0,ss,xloc,sig,shocks,order,compl,cond_shocks_id)

% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

narginchk(7,9)

n=nargin;
if n==7||isempty(compl)
    compl=@(x)1;
end
if n==7
    cond_shocks_id=[];
end
[nx,kplus1]=size(shocks);
if isempty(cond_shocks_id)
    cond_shocks_id=1:nx;
end

myzero=-sqrt(eps);
% badvector=@(x)any(compl(x)<myzero); % recreating this slows things down...
debug=false;
if debug
    fsolve_options=struct('Display','iter','TolFun',1e-12,'TolX',1e-12);
else
    fsolve_options=struct('Display','none','TolFun',1e-12,'TolX',1e-12);
end

y1=utils.forecast.one_step_engine(T,y0,ss,xloc,sig,shocks,order);

iter=0;
retcode=0;
% is_active_shock=false(1,kplus1);
if any(compl(y1.y)<myzero) % <---badvector(y1.y)
    nrows=size(compl(y1.y),1);
    if islogical(cond_shocks_id)
        cond_shocks_id=find(cond_shocks_id);
    end
    n_cond_shocks=numel(cond_shocks_id);
    if n_cond_shocks~=nrows
        error('Singularity in the problem of matching restrictions')
    end
    y1k=y0(:,ones(1,kplus1+1));
    done=false;
    while ~done
        iter=iter+1;
        % find the shocks that make the violations go away and the
        % constraints bind exactly
        %------------------------------------------------------------------
        shocks0=shocks;
        shocks0(cond_shocks_id,1:iter)=nan;
        e_id=isnan(shocks0);
        ee0=zeros(n_cond_shocks*iter,1);
        [ee,fval,exitflag]=fsolve(@multi_complementarity,ee0,fsolve_options);
        exitflag=utils.optim.exitflag(exitflag,ee,max(abs(fval(:))));
        if exitflag~=1
            retcode=701;
            is_active_shock=any(abs(shocks)>sqrt(eps),1); % is_active_shock(1:iter)=true;
            return
        end
        shocks0(e_id)=ee;
        % given those shocks, check that all future steps do not violate
        % the constraints
        %------------------------------------------------------------------
        test_passed=true;
        for icol=2:kplus1+1
            shocks_i=[shocks0(:,icol-1:end),zeros(nx,icol-2)];
            y1k(icol)=utils.forecast.one_step_engine(T,y1k(icol-1),ss,xloc,...
                sig,shocks_i,order);
            test_passed=~any(compl(y1k(icol).y)<myzero);% test_passed=~badvector(y1k(icol).y);
            if ~test_passed
                break
            end
            % here we could/should stop whenever there is a lift-up, since
            % there could be a chance that the step after the last one
            % still violates the constraints. But then this is true only if
            % the model is linear. If the model is nonlinear, it is
            % certainly a good idea to go to the end.
        end
        if iter==kplus1
            % the last guys is also constrained. Do one more step to see if
            % there is a lift up. If no lift up we are not out of the bush
            %--------------------------------------------------------------
            y1_last_uncond=utils.forecast.one_step_engine(T,y1k(end),ss,...
                xloc,sig,zeros(size(shocks0)),order);
            if ~all(compl(y1_last_uncond.y)>myzero)
                retcode=701;
                is_active_shock=any(abs(shocks)>sqrt(eps),1); % is_active_shock(1:iter)=true;
                return
            end
        end
        % if the test is passed we are done
        %----------------------------------
        if test_passed
            y1=y1k(2);
            done=true;
        end
    end
end
is_active_shock=any(abs(shocks)>sqrt(eps),1); % is_active_shock(1:iter)=true;

    function viol=multi_complementarity(ee)
        shocks0(e_id)=ee;
        viol=zeros(nrows,iter);
        for icol_=2:iter+1
            shocks_i=[shocks0(:,icol_-1:end),zeros(nx,icol_-2)];
            y1k(icol_)=utils.forecast.one_step_engine(T,y1k(icol_-1),ss,...
                xloc,sig,shocks_i,order);
            viol(:,icol_-1)=compl(y1k(icol_).y);
        end
        viol=viol(:);
    end

end