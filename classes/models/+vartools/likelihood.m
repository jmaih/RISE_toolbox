function varargout=likelihood(M,mapping,y,x,is_time_varying_trans_prob,...
    markov_chains)

% solve the model (more precisely the time-varying transition matrix) with
% the first observation 
first_obs=vartools.xpand_panel(y(:,1,:));

[sol,retcode]=vartools.solve(mapping,M,first_obs);

% number of solutions = number of parameter vectors = number of matrices
%-----------------------------------------------------------------------

nout=nargout;

outCell=cell(1,nout);

types=cell(1,nout);

ns=numel(sol);

for isol=1:ns
    
    [outCell{1:nargout}]=main_engine(sol(isol),M(:,:,isol),mapping,...
        y,x,is_time_varying_trans_prob,markov_chains,retcode(isol));
        
    if ns==1||isol==1
        
        varargout=outCell;
        
    end
    
    if ns>1
        
        restore()
        
    end
    
end

    function restore()
        
        for ii=1:nout
            
            outi=outCell{ii};
            
            if isol==1
                
                if isa(outi,'double')
                    
                    if isscalar(outi)
                        
                        types{ii}='d';
                        
                        varargout{ii}=outi*ones(1,ns);
                        
                    else
                        
                        types{ii}='v';
                        
                        varargout{ii}=outi(:,ones(1,ns));
                        
                    end
                    
                elseif isstruct(outi)
                    
                    types{ii}='s';
                    
                    varargout{ii}=drill(varargout{ii},outi,ns,isol);
                end
                
            else
                
                switch types{ii}
                    
                    case 'd'
                        
                        varargout{ii}(isol)=outi;
                        
                    case 'v'
                        
                        varargout{ii}(:,isol)=outi;
                        
                    case 's'
                        
                        varargout{ii}=drill(varargout{ii},outi,ns,isol);
                        
                end
                
            end
            
        end
        
    end

end
    
function y=drill(y,x,ns,sol_id)

if isstruct(x)
    
    fx=fieldnames(x);
    
    for ix=1:numel(fx)
        
        y.(fx{ix})=drill(y.(fx{ix}),x.(fx{ix}),ns,sol_id);
        
    end
    
else
    
    if sol_id==1
        
        y=x(:,:,ones(ns,1));
        
    else
        
        y(:,:,sol_id)=x;
        
    end       
    
end

end

function [LogLik,Incr,retcode,filtering]=main_engine(sol,M,mapping,...
    y,x,is_time_varying_trans_prob,markov_chains,retcode)

big_penalty=uminus(1e+8);

if retcode
    
    LogLik=big_penalty;
    
    Incr=[];
    
    filtering=[];
    
    return
    
end

h=mapping.nregimes;

B=sol.B;

SIG=sol.S;

Q0=sol.Q.Q;

y=y(:,:); % stretch for panel

x=x(:,:); % stretch for panel

[n,T]=size(y);

In=eye(n);

Q=Q0(:,:,ones(T+1,1));

% before seeing the data...
Q(:,:,1)=1/h;

% after seeing the data
if is_time_varying_trans_prob
    
    for tau=1:T
        
        fQ=mapping.solution.transition_matrix(M(:,1),y(:,tau));
        
        Q(:,:,tau+1)=fQ.Q;
        
    end
    
end

[Incr,ut,prob_update,prob_filt]=msvar_likelihood(B,SIG);
% sample to be used in the computation of the likelihood
%-------------------------------------------------------
start=1;

last=T;

% form the likelihood through summing the period-by-period densities.
% Note that here we return the value of the log-likelihood (for
% maximization). The negative will be taken inside estimate.m for
% minimization
%-------------------------------------------------------------------
LogLik=sum(Incr(start:last));

if retcode
    
    LogLik=big_penalty;
    
end

% negative of likelihood: will be taken inside estimate
%------------------------------------------------------
% % % %     LogLik=-LogLik;

if nargout<4
    
    return
    
end

regime_names=mapping.regimes(2:end,1).';

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
nshocks=size(ut{1},1);

exo_names=cell(1,nshocks);

for ishock=1:nshocks
    
    exo_names{ishock}=sprintf('shock_%0.0f',ishock);
    
end

filtering.smoothed_shocks=store_item(exo_names,ut);

Expected_smoothed_shocks=utils.filtering.expectation(prob_smooth,ut,true);

filtering.Expected_smoothed_shocks=store_item(exo_names,...
    {Expected_smoothed_shocks});

    function regime_filtration()
        
        filtering.filtered_regime_probabilities=store_prob(prob_filt);
        
        filtering.updated_regime_probabilities=store_prob(prob_update);
        
        prob_smooth=msvar_smoother();
        
        filtering.smoothed_regime_probabilities=store_prob(prob_smooth);
        
        function out=store_prob(item)
            
            out=struct();
            
            for st=1:h
                
                out.(regime_names{st})=item(st,:);
                
            end
            
        end
        
    end

    function state_filtration()
        
        state_names=markov_chains(1).state_names;
        
        for ic=2:numel(markov_chains)
            
            state_names=[state_names,markov_chains(ic).state_names]; %#ok<AGROW>
            
        end
        
        nstates=numel(state_names);
        
        nobs=size(prob_update,2);
        
        filt=zeros(nstates,nobs+1);
        
        updated=zeros(nstates,nobs);
        
        smooth=zeros(nstates,nobs);
        
        iter=0;
        
        regimes=cell2mat(mapping.regimes(2:end,2:end));
        
        nstates=numel(state_names);
        
        check_states=cell(1,nstates);
        
        for ic=1:numel(markov_chains)
            
            chain_name=markov_chains(ic).name;
            
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
        
        filtering.filtered_state_probabilities=...
            store_state_probabilities(filt);
        
        filtering.updated_state_probabilities=...
            store_state_probabilities(updated);
        
        filtering.smoothed_state_probabilities=...
            store_state_probabilities(smooth);
        
        function out=store_state_probabilities(datta)
            
            for ivar=1:nstates
                
                out.(check_states{ivar})=datta(ivar,:);
                
            end
            
        end
        
    end

    function out=store_item(names,item)
        
        out=struct();
        
        t0=size(item{1},2);
        
        hbar=min(h,numel(item));
        
        datta=nan(hbar,t0);
        
        for ix=1:numel(names)
            
            for st=1:hbar
                
                datta(st,:)=item{st}(ix,:);
                
            end
            
            if hbar==h
                
                out.(names{ix})=datta;
                
            else
                
                out.(names{ix})=datta;
                
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
        
        detSIG=nan(1,h);
        
        for st=h:-1:1 % reverse the order of computation
            
            [~,p] = chol(SIG(:,:,st));
            
            if p~=0
                
                retcode=305;
                
                return
                
            end
            
            ds=det(SIG(:,:,st));
            
            iSIG(:,:,st)=SIG(:,:,st)\In;
            
            detSIG(st)=((2*pi)^n*ds)^(-0.5);
            
            % compute the density conditional on the regime (st)
            %---------------------------------------------------
            ut{st}=y-B(:,:,st)*x;
            
        end
        
        exponent=zeros(1,h);
        
        for t=1:T
            % aggregate/integrate the densities
            %----------------------------------
            for st=h:-1:1 % reverse the order of computation
                
                ust=ut{st}(:,t);
                
                exponent(st)=-0.5*ust'*iSIG(:,:,st)*ust;
                
            end
            
            largest_exponent=max(exponent);
            
            for st=h:-1:1 % reverse the order of computation
                
                fst=detSIG(st)*exp(exponent(st)-largest_exponent);
                
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
            %             update_transition();
            
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
        
        %         function update_transition()
        %
        %             Qnext=Q(:,:,t);
        %
        %             Q(:,:,t+1)=Qnext;
        %
        %         end
        
    end

end