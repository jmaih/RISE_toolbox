function [T,R,steady_state,new_order,state_vars_location,log_vars]=...
load_solution(obj,type,do_bvar_dsge)
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

% when the type is ov, R is empty while T is of size solve_order x
% regimes_number. Basically, T is returned as Tz, Tzz, Tzzz, etc. where the
% rows that were originally in alphabetical order are put in order of
% simulation/processing
% when the type is ov, R corresponds to the shocks and T corresponds to the
% autoregressive elements and the matrices in T are square

if isempty(obj)
    
    if nargout>1
        
        error([mfilename,':: number of output arguments cannot exceed 1 when the object is empty'])
    
    end
    
    T=cell(0,4);
    
    return
    
end

if nargin<3
    
    do_bvar_dsge=false;
    
end

if ~any(strcmp(type,{'ov','iov'}))
    
    error('type should be ov (order var) or iov (inv order var)')
    
end

is_alphabetical_order=strcmp(type,'iov');

if do_bvar_dsge
    
    if ~obj.is_dsge_var_model && obj.options.dsgevar_var_regime
        
        error('requesting the VAR solution without conditions being met')
        
    end
    
    is_alphabetical_order=true;
    
end

regimes_number=obj.markov_chains.regimes_number;

order=obj.options.solve_order;

T=cell(order,regimes_number);

R=cell(1,regimes_number);

ov=obj.order_var;

steady_state=obj.solution.ss;

new_order=ov;

state_vars_location=obj.locations.after_solve.t.pb;

log_vars=obj.log_vars;
        
% update T and steady state
%--------------------------
order_var_solution();

% now get the inv_order_var solution if necessary
%------------------------------------------------
if is_alphabetical_order
    
    iov=obj.inv_order_var;
    
    % only for order 1
    inv_order_var_solution()
    
    new_order=1:numel(new_order);
    
    state_vars_location=[];
    
end
    
if do_bvar_dsge
    
    k_minus_constant=bvar_dsge_solution();
    
    new_order=1:k_minus_constant;
    
end

    function k_minus_constant=bvar_dsge_solution()
        
        T_bvardsge=obj.dsge_var.T;
        
        n=obj.dsge_var.n;
        
        p=obj.dsge_var.p;
        
        constant=obj.dsge_var.constant;
        
        k_minus_constant=n*p;
        
        PHI=obj.dsge_var.posterior.PHI;
        
        SIG=obj.dsge_var.posterior.SIG;
        
        if obj.options.dsgevar_inner_param_uncertainty
            
            lambda=obj.dsge_var.lambda;
            
            S=(1+lambda)*T_bvardsge*SIG;
            
            df=obj.dsge_var.posterior.inverse_wishart.df(1);
            
            [~,~,~,rndfn]=distributions.inv_wishart();
            
            SIG_draw=rndfn(S,df);
            
            COV_PHI=kron(SIG_draw,obj.dsge_var.posterior.ZZi);
            
            CS=transpose(chol(COV_PHI));
            
            is_explosive=true;
            
            while is_explosive
                
                PHI_draw=PHI(:)+CS*randn(size(CS,1),1);
                
                PHI_draw=reshape(PHI_draw,size(PHI));
                
                T{1}=companion_form(PHI_draw(constant+1:end,:)');
                
                is_explosive=max(eig(T{1}))>=obj.options.stability_criterion;
                
            end
            
        else
            
            PHI_draw=PHI;
            
            SIG_draw=SIG;
            
            T{1}=companion_form(PHI_draw(constant+1:end,:)');
            
        end
        
        % compute the mean to be used in various simulations?
        if constant
            
            c=[PHI_draw(1,:)';zeros(n*(p-1),1)];
            
            mu=(eye(k_minus_constant)-T{1})\c;
            
        else
            mu=zeros(k_minus_constant,1);
            
        end
        
        steady_state{1}=mu;
        
        % Get rotation
        stochastic=~obj.exogenous.is_observed;
        
        nx=sum(stochastic);
        
        Atheta = R{1}(obj.observables.state_id,stochastic); %
        
        [OMEGAstar,SIGMAtr] = qr(Atheta');  %#ok<ASGLU>
        
        % identification
        %----------------
        SIGtrOMEGA = transpose(chol(SIG_draw))*OMEGAstar';
        
        R{1}=[SIGtrOMEGA;zeros(n*(p-1),nx)];% R{1}=SIGtrOMEGA;% 
        
        % splice and destroy
        %-------------------
        T{1}=[T{1},zeros(k_minus_constant,1),R{1}]; %[T,sig,R]
        
        function c=companion_form(x)
            
            c=[x
            eye(n*(p-1)),zeros(n*(p-1),n)];
        
        end
        
    end

    function inv_order_var_solution()
        
        endo_nbr_=obj.endogenous.number;
        
        exo_nbr_=sum(obj.exogenous.number);
        
        z_pb=obj.locations.after_solve.z.pb;
        
        t_pb=state_vars_location;
        
        e_0=obj.locations.after_solve.z.e_0;
        
        Tz=R;
        
        tmp=zeros(endo_nbr_);
        
        for isol=1:regimes_number
            % this solution comes from order_var below and so we re-order
            % both rows and columns
            tmp(:,t_pb)=T{1,isol}(:,z_pb);
            
            % separate autoregressive part from shocks
            %-----------------------------------------
            Tz{isol}=tmp(iov,iov);
            
            R{isol}=T{1,isol}(iov,e_0(1):end);
            
            if regimes_number>1
                
                if isol==1
                    
                    npges=size(R{isol},2)/exo_nbr_;
                    
                end
                
                R{isol}=reshape(full(R{isol}),[endo_nbr_,exo_nbr_,npges]);
                
            end
            
        end
        
        T=Tz;
        
        if order>1
            
            warning([mfilename,':: Only first order will be used'])
            
        end
        
        % re-order the log vars
        %----------------------
        log_vars=log_vars(iov);
    end

    function order_var_solution()
        
        zzz=repmat('z',1,order);
        
        is_dev=obj.options.simul_bgp_deviation;
        
        target=numel(state_vars_location)+1;
        
        for io=1:order
            
            for ireg=1:regimes_number
                
                if io==1
                    % Take log of log-var variables to be consistent with
                    % the solution
                    %-----------------------------------------------------
                    steady_state{ireg}(log_vars)=...
                        log(steady_state{ireg}(log_vars));
                    
                    if ~is_alphabetical_order
                        % do this only if the alphabetical order is not
                        % needed 
                        steady_state{ireg}=steady_state{ireg}(ov);
                        
                    end
                    
                end
                % re-order the rows
                %------------------
                T{io,ireg}=obj.solution.(['T',zzz(1:io)]){ireg}(ov,:);
                
                if is_dev && io==1
                    
                    T{io,ireg}(:,target)=real(T{io,ireg}(:,target));
                    
                end
                
            end
            
        end
        
        % Now also re-order the log vars
        %-------------------------------
        log_vars=log_vars(ov);
        
    end

end