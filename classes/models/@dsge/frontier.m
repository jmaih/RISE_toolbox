function f=frontier(obj,lambda_name,lambda_vals,simul)
% frontier -- computes standard devations of the model for a grid over a
% given parameter
%
% Syntax
% -------
% ::
%
%   f=frontier(obj,lambda_name,lambda_vals)
%
% Inputs
% -------
%
% - **obj** [rise|dsge]: model object
%
% - **lambda_name** [char]: name of the parameter to vary
%
% - **lambda_vals** [vector]: 1 x 2 or 1 x N vector of values for the
% parameter. When N=2, a grid of 50 points is constructed between the two
% values. When N>2, lambda_vals is the grid.
%
% - **simul** [true|{false}]: use simulation instead of theoretical
% moments.
%
% Outputs
% --------
%
% - **f** [struct]: standard deviations of all variables in the model for
% each value of lambda_name
%
% More About
% ------------
%
% - improvements to consider are how to deal with regime switches or
% nonlinear models in general. One solution is to use simulation
%
% Examples
% ---------
%
% See also:

if isempty(obj)
    f=struct();
    return
end

if nargin<4
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

nworkers=get_number_of_workers();
parfor(ival=1:nvals,nworkers)
    [V,retcode]=objective(lambda_vals(ival));
    if ~retcode
        sd(:,ival)=sqrt(diag(V));
    end
end

f=struct();
for iname=1:n
    f.(names{iname})=sd(iname,:);
end
f.lambda__=lambda_vals;
f.ngrid__=nvals;
f.simul_periods__=0;
if simul
    f.simul_periods__=obj.options.simul_periods;
end

    function n=get_number_of_workers()
        v=ver;
        Names={v.Name};
        v=v(strcmp(Names,'MATLAB'));
        v=str2double(v.Version);
        n=0;
        if v>8.1
            pool = gcp('nocreate');
            if ~isempty(pool)
                n=pool.NumWorkers;
            end
        else
            n=matlabpool('size'); %#ok<DPOOL>
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
            % use the same seed at the beginning of each simulation
            %-------------------------------------------------------
            [db,~,retcode]=simulate(obj,...
                'parameters',{lambda_name,val});
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