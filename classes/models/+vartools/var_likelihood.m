function [LogLik,Incr,retcode,obj]=var_likelihood(params,obj)%

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
        filtering.filtered_probabilities=store_prob(prob_filt);
        filtering.updated_probabilities=store_prob(prob_update);
        filtering.smoothed_probabilities=store_prob(msvar_smoother());
        exo_names=obj.exogenous.name(~obj.exogenous.is_observed);
        filtering.smoothed_shocks=store_item(exo_names,ut);
        obj=set(obj,'filters',filtering);
    end
    
    if obj.options.debug
        disp(LogLik)
    end
end

    function out=store_prob(item)
        out=struct();
        for st=1:h
            out.(regime_names{st})=ts(start_date,item(st,:)');
        end
    end
    function out=store_item(names,item)
        out=struct();
        t0=size(item{1},2);
        datta=nan(t0,h);
        for ix=1:numel(names)
            for st=1:h
                datta(:,st)=item{st}(ix,:)';
            end
            out.(names{ix})=ts(start_date,datta,regime_names);
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
        %         fst=cell(1,h);
        Incr=nan(T,1);
        pf=nan(h,1);
        prob_udate=zeros(h,T);
        prob_filt=zeros(h,T+1);
        [prob_filt(:,1),retcode]=initial_markov_distribution(Q(:,:,1),true);
        if retcode
            return
        end
        
        iSIG=SIG;
        detSIG=SIG;
        for st=1:h
            [~,p] = chol(SIG{st});
            if p~=0
                retcode=305;
                return
            end
            ds=det(SIG{st});
            iSIG{st}=SIG{st}\eye(n);
            detSIG{st}=((2*pi)^n*ds)^(-0.5);
        end
        % compute the density conditional on the regime (st)
        %---------------------------------------------------
        for st=1:h
            ut{st}=y-B{st}*x;
        end
        
        for t=1:T
            % aggregate/integrate the densities
            %----------------------------------
            for st=1:h
                ust=ut{st}(:,t);
                fst=detSIG{st}*exp(-0.5*ust'*iSIG{st}*ust);
                pf(st)=prob_filt(st,t)*fst;
            end
            lik_t=sum(pf);
            if ~isreal(lik_t)||lik_t<=0
                retcode=306; % unlikely parameter vector
                return
            end
            % compute the updated probabilities
            %----------------------------------
            prob_udate(:,t)=pf/lik_t;
            
            % add an increment to the likelihood
            %-----------------------------------
            Incr(t)=log(lik_t);
            
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