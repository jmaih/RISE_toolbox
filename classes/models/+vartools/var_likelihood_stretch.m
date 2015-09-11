function [LogLik,Incr,retcode,obj]=var_likelihood_stretch(params,obj)%
% VAR_LIKELIHOOD_STRETCH - computes the likelihood of rfvar and svar
% models.
%
% Syntax
% -------
% ::
%   [LogLik,Incr,retcode,obj]=VAR_LIKELIHOOD_STRETCH(params,obj)
%
% Inputs
% -------
%
% - **params** [empty|vector]: parameter vector under estimation
%
% - **obj** [rfvar|svar]: model object
%
% Outputs
% --------
%
% - **LogLik** [numeric]: log of the likelihood
%
% - **Incr** [vector]: contributions of each period/iteration to the
% likelihood
%
% - **retcode** [numeric]: 0 if no problem
%
% - **obj** [svar|rfvar]: model object possibly updated with the solution
% of the model and the filters (shocks and probabilities).
%
% More About
% ------------
%
% - For numerical reasons, this function STRETCHES the exponential function
% in order to avoid early truncations for low levels of the likelihood.
%
% Examples
% ---------
%
% See also: VARTOOLS/VAR_LIKELIHOOD, VARTOOLS/VAR_LIKELIHOOD_DIRECT

if isempty(obj)
    LogLik=struct();
    return
end

if nargin~=2
    error([mfilename,':: Number of arguments must be 2'])
end
is_endogenous_switching=false;
obj=assign_estimates(obj,params);

[obj,retcode]=solve(obj);
if retcode
    LogLik=obj.options.estim_penalty;
    Incr=[];
    return
end
h=obj.markov_chains.regimes_number;
nlags=obj.nlags;
[B,SIG]=vartools.resolve(obj.solution,nlags,h);
Q0=obj.solution.transition_matrices.Q;

[y,x,n,T]=vartools.set_y_and_x(obj.data.y,obj.data.x,nlags,obj.constant);

Q=nan(h,h,T+1);
Q(:,:,1)=Q0;

[Incr,ut,prob_update,prob_filt]=msvar_likelihood(B,SIG);
% sample to be used in the computation of the likelihood
%-------------------------------------------------------
start=1;
last=T;
if retcode
    % negative of likelihood will be taken inside estimate
    %-----------------------------------------------------
    LogLik=uminus(obj.options.estim_penalty);
    
    if obj.options.debug
        bmle=y/x;
        smle=y-bmle*x;
        smle=(smle*smle')/T;
        utils.error.decipher(retcode)
        [Incr]=msvar_likelihood({bmle},{smle});
        fprintf(1,'log-likelihood (mle) %0.0f\n',sum(Incr(start:last)));
    end
else
    
    % form the likelihood through summing the period-by-period densities.
    % Note that here we return the value of the log-likelihood (for
    % maximization). The negative will be taken inside estimate.m for
    % minimization
    %-------------------------------------------------------------------
    LogLik=sum(Incr(start:last));
    
    % negative of likelihood: will be taken inside estimate
    %------------------------------------------------------
    % % % %     LogLik=-LogLik;
    
    if ~obj.estimation_under_way
        start_date=obs2date(obj.options.estim_start_date,nlags+1);
        regime_names=obj.markov_chains.regime_names;
        filtering=struct();
        prob_smooth=[];
        % regimes filtration
        %---------------------
        regime_filtration()
        % state filtration
        %------------------
        state_filtration()
        
        % shocks and storage
        %--------------------
        exo_names=obj.exogenous.name(~obj.exogenous.is_observed);
        filtering.smoothed_shocks=store_item(exo_names,ut);
        Expected_smoothed_shocks=utils.filtering.expectation(prob_smooth,ut,true);
        filtering.Expected_smoothed_shocks=store_item(exo_names,...
            {Expected_smoothed_shocks});
        obj=set(obj,'filters',filtering);
    end
    
    if obj.options.debug
        disp(LogLik)
    end
end

    function regime_filtration()
        filtering.filtered_regime_probabilities=store_prob(prob_filt);
        filtering.updated_regime_probabilities=store_prob(prob_update);
        prob_smooth=msvar_smoother();
        filtering.smoothed_regime_probabilities=store_prob(prob_smooth);
        function out=store_prob(item)
            out=struct();
            for st=1:h
                out.(regime_names{st})=ts(start_date,item(st,:)');
            end
        end
    end

    function state_filtration()
        nstates=numel(obj.markov_chains.state_names);
        nobs=size(prob_update,2);
        filt=zeros(nstates,nobs+1);
        updated=zeros(nstates,nobs);
        smooth=zeros(nstates,nobs);
        iter=0;
        regimes=cell2mat(obj.markov_chains.regimes(2:end,2:end));
        state_names=obj.markov_chains.state_names;
        nstates=numel(state_names);
        check_states=cell(1,nstates);
        for ic=1:obj.markov_chains.chains_number
            chain_name=obj.markov_chains.chain_names{ic};
            max_states=max(regimes(:,ic));
            for istate=1:max_states
                iter=iter+1;
                check_states{iter}=sprintf('%s_%0.0d',chain_name,istate);
                this_regimes=find(regimes(:,ic)==istate);
                for ireg=1:numel(this_regimes)
                    filt(iter,:)=filt(iter,:)+prob_filt(this_regimes(ireg),:);
                    updated(iter,:)=updated(iter,:)+prob_update(this_regimes(ireg),:);
                    smooth(iter,:)=smooth(iter,:)+prob_smooth(this_regimes(ireg),:);
                end
            end
        end
        
%         disp('Any differences?'),disp(check_states),disp(state_names)
%         keyboard
        
        start_date=get(filtering.filtered_regime_probabilities.regime_1,'start');
        filtering.filtered_state_probabilities=...
            store_state_probabilities(filt.');
        filtering.updated_state_probabilities=...
            store_state_probabilities(updated.');
        filtering.smoothed_state_probabilities=...
            store_state_probabilities(smooth.');
        function out=store_state_probabilities(datta)
            for ivar=1:nstates
                if ivar==1
                    out.(check_states{ivar})=ts(start_date,datta(:,ivar));
                else
                    out.(check_states{ivar})=...
                        reset_data(out.(check_states{ivar-1}),datta(:,ivar));
                end
            end
        end
    end

    function out=store_item(names,item)
        out=struct();
        t0=size(item{1},2);
        hbar=min(h,numel(item));
        datta=nan(t0,hbar);
        for ix=1:numel(names)
            for st=1:hbar
                datta(:,st)=item{st}(ix,:)';
            end
            if hbar==h
                out.(names{ix})=ts(start_date,datta,regime_names);
            else
                out.(names{ix})=ts(start_date,datta);
            end
        end
    end

    function prob_smooth=msvar_smoother()
        prob_smooth=prob_update;
        for t=T-1:-1:1
            for st=1:h
                prob_smooth(st,t)=prob_update(st,t)*...
                    sum(prob_smooth(:,t+1)'.*Q(st,:,t+1)./prob_filt(:,t+1)');
                % here I need to check whether
                % - prob_filt should be t or t+1
                % - Q should be t or t+1
            end
        end
    end

    function [Incr,ut,prob_udate,prob_filt]=msvar_likelihood(B,SIG)
        ut=cell(1,h);
        Incr=nan(T,1);
        fake_pf=zeros(h,1);
        prob_udate=zeros(h,T);
        prob_filt=zeros(h,T+1);
        [prob_filt(:,1),retcode]=initial_markov_distribution(Q(:,:,1),true);
        if retcode
            return
        end
        
        iSIG=SIG;
        detSIG=SIG;
        for st=h:-1:1 % reverse the order of computation
            [~,p] = chol(SIG{st});
            if p~=0
                retcode=305;
                return
            end
            ds=det(SIG{st});
            iSIG{st}=SIG{st}\eye(n);
            detSIG{st}=((2*pi)^n*ds)^(-0.5);
            
            % compute the density conditional on the regime (st)
            %---------------------------------------------------
            ut{st}=y-B{st}*x;
        end
        
        exponent=zeros(1,h);
        for t=1:T
            % aggregate/integrate the densities
            %----------------------------------
            for st=h:-1:1 % reverse the order of computation
                ust=ut{st}(:,t);
                exponent(st)=-0.5*ust'*iSIG{st}*ust;
            end
            largest_exponent=max(exponent);
            for st=h:-1:1 % reverse the order of computation
                fst=detSIG{st}*exp(exponent(st)-largest_exponent);
                fake_pf(st)=prob_filt(st,t)*fst;
            end
            fake_lik_t=sum(fake_pf);
            % lik_t=fake_lik_t*exp(largest_exponent) sometimes too large
            if ~isreal(fake_lik_t)||fake_lik_t<=0
                retcode=306; % unlikely parameter vector
                return
            end
            % compute the updated probabilities
            %----------------------------------
            prob_udate(:,t)=fake_pf/fake_lik_t;
            % =pf/lik_t=(exp(largest_exponent)*fake_pf)/(exp(largest_exponent)*fake_lik_t);
            
            % add an increment to the likelihood
            %-----------------------------------
            Incr(t)=largest_exponent+log(fake_lik_t);
            
            % update transition if necessary
            %-------------------------------
            update_transition();
            if retcode
                return
            end
            
            % compute the filtered probabilities to be used for next
            % period's aggregation of the conditional likelihoods
            %-------------------------------------------------------
            for st=1:h
                for slag=1:h
                    prob_filt(st,t+1)=prob_filt(st,t+1)+Q(slag,st,t+1)*prob_udate(slag,t);
                end
            end
        end
        function update_transition()
            if is_endogenous_switching
                error('please remind junior.maih@gmail.com to implement endogenous switching for vars')
                Qnext=something(y(:,t),params);
                if ~all(sum(Qnext,2)==1)
                    retcode=3; % problem in transition matrix
                end
            else
                Qnext=Q(:,:,t);
            end
                Q(:,:,t+1)=Qnext;
        end
    end
end