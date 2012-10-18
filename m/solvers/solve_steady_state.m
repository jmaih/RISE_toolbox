function [ys,retcode]=solve_steady_state(ys0,x_ss,param,resid_func,static_func,is_linear_model,optim_opt)

last_item=numel(ys0);
switch is_linear_model
    case 1
        ys0=0*ys0;
        % compute the constant
        const=resid_func(ys0,static_func,x_ss,param);
        if any(const)
            % partial derivatives
            PD=zeros(last_item);
            for jj=1:last_item
                yj=ys0;yj(jj)=1;
                PD(:,jj)=resid_func(yj,static_func,x_ss,param)-const;
            end
            ys=-pinv(PD)*const; % <-- ys=-PD\const;
            % The generalized inverse works better when the system is
            % singular as it is the case when the model is solved with
            % unit roots. But it might be slower and so I don't know
            % whether I should put a switch...
        else
            ys=ys0;
        end
        residuals=resid_func(ys,static_func,x_ss,param);
        exitflag=max(abs(residuals))<=optim_opt.TolFun;
    otherwise
        % Here is why we need good initial values. You cannot, say start at
        % zero for a variable in logs...
        residuals=resid_func(ys0,static_func,x_ss,param);
        if isnan(residuals)
            exitflag=-4;
            ys=nan;
        elseif max(abs(residuals))<=optim_opt.TolFun
            ys=ys0;
            exitflag=1;
        else
            optim_opt.Display='none';
            [ys,junk,exitflag]=fsolve(resid_func,ys0,optim_opt,static_func,x_ss,param);
            if exitflag~=1
                residuals=resid_func(ys,static_func,x_ss,param);
                if ~all(isnan(residuals))&& max(abs(residuals))<=optim_opt.TolFun
                    exitflag=1;
                end
            end
        end
        %         [ys,RESNORM,RESIDUAL,exitflag]=lsqnonlin(resid_func,ys0,[],[],optim_opt,static_func,x_ss,param,def);
end
if ~any(isnan(ys)) && exitflag==1
    % maybe I should not do the following? especially here?
    ys(abs(ys)<1e-12)=0;
    retcode=0;
else
    retcode=1; % steady state cannot be solved
end

