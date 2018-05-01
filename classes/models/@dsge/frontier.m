function f=frontier(obj,lambda_name,lambda_vals,simul,seed)
% frontier -- computes standard devations of the model for a grid over a
% given parameter
%
% ::
%
%
%   f=frontier(obj,lambda_name,lambda_vals)
%
% Args:
%
%    - **obj** [rise|dsge]: model object
%
%    - **lambda_name** [char]: name of the parameter to vary
%
%    - **lambda_vals** [vector]: 1 x 2 or 1 x N vector of values for the
%    parameter. When N=2, a grid of 50 points is constructed between the two
%    values. When N>2, lambda_vals is the grid.
%
%    - **simul** [true|{false}]: use simulation instead of theoretical
%    moments.
%
%    - **seed** [numeric|{1971}]: see for simulations
%
% Returns:
%    :
%
%    - **f** [struct]: standard deviations of all variables in the model for
%    each value of lambda_name and some further information about the
%    simulation process under a substructure with name **stats__**. The fields
%    of the sub-structure are:
%      - **lambda** [vector]: discretized values of lambda
%      - **ngrid** [scalar]: number of grid points
%      - **simul_periods** [integer]: number of simulations. If strictly
%      positive, then simulation is used for computing the moments of the
%      process
%      - **retcode** [vector]: information on how successful each run in the
%      grid is. if retcode=0, there is no problem. If retcode different from
%      zero, then the information can be retrieved by running
%      decipher(retcode)
%
% Note:
%
%    - improvements to consider are how to deal with regime switches or
%    nonlinear models in general. One solution is to use simulation
%
% Example:
%
%    See also:

if isempty(obj)
    
    f=cell(0,4);
    
else
    
    if nargin<5
        
        seed=[];
        
        if nargin<4
            
            simul=[];
            
        end
        
    end
    
    if isempty(seed)
        
        seed=1971;
        
    end
    
    
    if isempty(simul)
        
        simul=false;
        
    end
    
    lambda_vals=sort(lambda_vals);
    
    if numel(lambda_vals)==2
        
        lambda_vals=linspace(lambda_vals(1),lambda_vals(2),50);
        
    end
    
    is_linear=obj.options.solve_order==1 && obj.markov_chains.regimes_number;
    
    if ~isempty(simul)
        
        if ~islogical(simul)
            
            error('simul must be empty or logical')
            
        end
        
    end
    
    simul=simul||~is_linear;
    
    nvals=numel(lambda_vals);
    
    obj=set(obj,'autocov_ar',0,...
        'simul_to_time_series',false);%,'lyapunov_algo','schur'
    
    n=obj.endogenous.number;
    
    names=obj.endogenous.name;
    
    good_locs=[];
    
    objective=@(x)variance_engine(x);
    
    sd=nan(n,nvals);
    
    nworkers=utils.parallel.get_number_of_workers();
    
    retcode=zeros(1,nvals);
    
    parfor(ival=1:nvals,nworkers)
        
        [V,retcode(ival)]=objective(lambda_vals(ival));
        
        if ~retcode(ival)
            
            sd(:,ival)=sqrt(diag(V));
            
        end
        
    end
    
    f=struct();
    
    for iname=1:n
        
        f.(names{iname})=sd(iname,:);
        
    end
    
    f.stats__=struct('lambda',lambda_vals,'ngrid',nvals,'simul_periods',0,...
        'retcode',retcode);
    
    if simul
        
        f.stats__.simul_periods__=obj.options.simul_periods;
        
    end
    
end

    function [V,retcode]=variance_engine(val)
        
        if simul
            
            [V,retcode]=simulated_variances();
            
        else
            
            [V,retcode]=theoretical_autocovariances(obj,...
                'parameters',{lambda_name,val});
            
        end
        
        function [V,retcode]=simulated_variances()
            
            % reset random number generator
            %------------------------------
            rng(seed)
            
            % use the same seed at the beginning of each simulation
            %-------------------------------------------------------
            [db,~,retcode]=simulate(obj,'parameters',{lambda_name,val});
            
            if retcode
                
                V=[];
                
            else
                % which one of the workers will fill out good_locs first?
                % All of them: good_locs seems to become a local variable
                % to each slave's workspace
                %--------------------------------------------------------
                if isempty(good_locs)
                    
                    good_locs=locate_variables(names,db{2}{2,2});
                    
                end
                
                V=cov(db{1}(:,good_locs));
                
            end
            
        end
        
    end

end