function db=randsample(obj,N,K,varargin)
% RANDSAMPLE -- generates random time series from a model
%
% ::
%
%
%   db=randsample(obj,N,K)
%
%   db=randsample(obj,N,K,varargin)
%
% Args:
%
%    - **obj** [rise|dsge|svar|rfvar]: model object
%
%    - **N** [integer]: length of the time series
%
%    - **K** [integer]: number of time series
%
%    - **varargin** []: additional options for the model object
%
% Returns:
%    :
%
%    - **db** [struct]: time series
%
% Note:
%
% Example:
%
%    See also:


% returns K samples of length N

% The solution for do not anticipate should be killed inside initial
% conditions

if isempty(obj)
    
    db=cell(0,4);
    
    return
    
end

[obj,retcode]=solve(obj,varargin{:});

if retcode
    
    error(['model cannot be solved:: ',decipher(retcode)])
    
end

% initial conditions
%-------------------
[y0,ss,T,xloc,sig,order,compl,Qfunc,is_log_var,endHist,k_future]=load_initial_conditions();

nx=sum(obj.exogenous.number);
            
shocks=zeros(nx,k_future+1,N,K);
            
y=[]; regimes=nan(K,N);

h=obj.markov_chains.regimes_number;
            
for k=1:K
    
    y00=y0;
    
    t=0;
    
    redo_Q=true;
    
    while t<N
        
        t=t+1;
        
        shocks_t=randn(nx,k_future+1);
        
        % compute transition matrix and switching probabilities
        %------------------------------------------------------
        if redo_Q
            
            rt=regimes(k,t);
            
            [Q,retcode]=Qfunc(y00.y);
            
            if retcode
                % continue on other simulations...
                break
                
            end
            
            state_list=1:h;
            
            % the state is not known
            if t==1
                % draw from initial distribution
                PAI=1/h*ones(1,h);
            else
                % draw conditional on yesterday's state
                PAI=Q(regimes(k,t-1),:);
            end
            
        end
        
        rt=draw_regime(rt);
        
        if isnan(rt)
            % feasible regime could not be found
            
            break
            
        end
        
        y1=utils.forecast.one_step_engine(T(:,rt),y00,ss{rt},...
            xloc,sig,shocks_t,order);
        
        if isfeasible(y1)
            
            redo_Q=true;
            
            store_output();
            
            y00=y1;
            
        else
            
            redo_Q=false;
            
            t=t-1;
            
            % regime fails and is terminated
            state_list(state_list==rt)=[];
            
            % reset the regime
            rt=nan;
            
        end
         
    end
   
end

% format output
%---------------
db=format_output();

    function db=format_output()
        % add initial conditions
        %------------------------
        y=cat(2,y0.y(:,:,ones(1,K)),y);
        y0cols=size(y0.y,2);
        new_N=size(y,2);
        
        shocks=cat(3,zeros(nx,k_future+1,y0cols,K),shocks);
        
        regimes=cat(2,nan(K,y0cols),regimes);
        
        % exponentiate before going any further: is_log_var is in the
        % order_var order
        %------------------------------------------------------------------
        if ~isempty(is_log_var) && any(is_log_var)
            
            y(is_log_var,:,:)=exp(y(is_log_var,:,:));
            
        end
        
        % put y in the correct order before storing
        %------------------------------------------
        y=re_order_output_rows(obj,y);
        
        % time x variable x replications
        %--------------------------------
        y=permute(y,[2,1,3]);
        regimes=permute(regimes,[2,3,1]);
        
        % time x variable x replications x horizon
        %-----------------------------------------
        shocks=permute(shocks,[3,1,4,2]);
        
        % find the corresponding states
        %-------------------------------
        states=regimes2states(obj,regimes);
        
        start_date=date2serial(endHist)-y0cols+1;
        
        db=struct();
        
        ysrs=cat(2,y,regimes,states);
        
        vnames=[get(obj,'endo_list'),'regime',...
            get(obj,'state_list')-'const_1'];
        
        replic_list=parser.create_state_list('replic',K);
        
        load_batch(ysrs,vnames,1)
        
        load_batch(shocks,get(obj,'exo_list'),k_future+1)
               
        function load_batch(data,vnames,npages)
            
            d=reshape(data(:,1,:,:),new_N,K,npages);
            % d=squeeze(data(:,1,:,:)); does not work well if K=1
            
            prototype=ts(start_date,d,replic_list,'',true);
            
            for ivar=1:numel(vnames)
                
                d=reshape(data(:,ivar,:,:),new_N,K,npages);
                
                db.(vnames{ivar})=reset_data(prototype,d);
                
            end
            
        end

    end

    function flag=isfeasible(x)
        
        flag=isempty(compl) || compl(x.y);
        
    end

    function reg=draw_regime(reg)
        
        if isnan(reg)
            
            cp=rebuild_cp();
            
            if isempty(cp)
                
                return
                
            end
            
            lucky=find(cp>rand,1,'first')-1;
            
            reg=state_list(lucky);
            
        end
        
        function cp=rebuild_cp()
            
            cp=[];
            
            if isempty(state_list)
                
                warning('I could not find a feasible path (all paths exhausted)')
                
                return
                
            end
            
            PAI00=PAI(state_list);
            
            PAI00=PAI00/sum(PAI00);
            
            if any(isnan(PAI00))
                
                warning('I could not find a feasible path (prob=0)')
                
                return
                
            end
            
            cp=cumsum(PAI00);
            
            cp=[0,cp(:).'];
            
        end
        
    end

    function store_output()
        
        if isempty(y)
            
            ny=numel(y1.y);
            
            y=zeros(ny,N,K);
            
        end
        
        regimes(k,t)=rt;
        
        y(:,t,k)=y1.y;
        
        shocks(:,:,t,k)=shocks_t;
        
    end

    function [y0,ss,T,xloc,sig,order,compl,Qfunc,is_log_var,endHist,k_future]=...
            load_initial_conditions()
        
        Initcond=set_simulation_initial_conditions(obj);
        
        y0=Initcond.y;
        
        if numel(y0)>1
            
            error('more than one initial conditions')
            
        end
        
        % use the steady state with possibly loglinear variables
        ss=Initcond.log_var_steady_state;
        
        T=Initcond.T;
        
        xloc=Initcond.state_vars_location;
        
        sig=Initcond.simul_sig;
        
        order=Initcond.simul_order;
        
        compl=Initcond.complementarity;
        
        Qfunc=Initcond.Qfunc;
        
        is_log_var=Initcond.is_log_var;
        
        endHist=Initcond.simul_history_end_date;
        
        k_future=Initcond.k_future;
        
    end

end

function states=regimes2states(obj,regimes_history)

[N,Wan,replic]=size(regimes_history);

if ~isequal(Wan,1)
    
    error('confusing ordering of regimes')
    
end

regimes_tables=obj.markov_chains.regimes;

markov_chains=regimes_tables(1,2:end);

nchains=numel(markov_chains);

reg_table=cell2mat(regimes_tables(2:end,2:end));

states=nan(N,nchains,replic);

for ip=1:replic
    % there is only one column and so the following should also
    % work max_reg=max(regimes_history(:,ip));

    max_reg=max(regimes_history(:,1,ip));
    
    for ireg=1:max_reg
        
        pos=regimes_history(:,ip)==ireg;
        
        n=sum(pos);
        
        states(pos,:,ip)=reg_table(ireg*ones(n,1),:);
        
    end
    
end

end