function [xparam,retcode,could_solve,behave,logpost,...
    LogLik,log_prior]=monte_carlo_filtering(obj,varargin)

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    xparam=struct('monte_carlo_nsim',2^12,'monte_carlo_behavioral_function',[]);
    return
end
nn=length(varargin);
nobj=numel(obj);
if nn
    if rem(nn,2)~=0
        error([mfilename,':: arguments must come in pairs'])
    end
    discard=[];
    for ii=1:.5*nn
        if strcmp(varargin{2*ii-1},'optimset')
            for jj=1:nobj
                obj(jj).options.optimset=utils.miscellaneous.setfield(obj(jj).options.optimset,varargin{2*ii});
            end
            discard=[2*ii-1,2*ii];
        end
    end
    varargin(discard)=[];
    if ~isempty(varargin)
        for jj=1:nobj
            obj(jj).options=utils.miscellaneous.setfield(obj(jj).options,varargin{:});
        end
    end
end

monte_carlo_behavioral_function=obj.options.monte_carlo_behavioral_function;
monte_carlo_nsim=obj.options.monte_carlo_nsim;

varargs={};
if ~isempty(monte_carlo_behavioral_function)
    if iscell(monte_carlo_behavioral_function)
        varargs=monte_carlo_behavioral_function(2:end);
        monte_carlo_behavioral_function=monte_carlo_behavioral_function{1};
    end
end

behavior=~isempty(monte_carlo_behavioral_function);
% pval_cutoff=0.1;

%% get the data for likelihood computation
smooth_flag=false;
if obj.observables.number(1)>0 % || isempty(obj.varobs)
    % load the data
    [obj,~,retcode]=load_data(obj);
    if retcode
        warning([mfilename,':: ',utils.error.decipher(retcode)]) %#ok<WNTAG>
    else
        smooth_flag=true; %~isempty(vertcat(obj.varobs.value));
%         if ~smooth_flag
%             warning([mfilename,':: data not available, likelihood won''t be computed']) %#ok<WNTAG>
%         end
    end
end

%% sample parameters: this is where I should call sobol, halton, latin_hypercubes, etc.
bounds=[vertcat(obj.estimation.priors.lower_bound),vertcat( obj.estimation.priors.upper_bound)];
npar=size(bounds,1);
xparam=quasi_monte_carlo.sobol(npar,monte_carlo_nsim,bounds(:,1),bounds(:,2)); % number of sub-intervals
% % xparam=bsxfun(@plus,bounds(:,1),bsxfun(@times,rand(npar,monte_carlo_nsim),bounds(:,2)-bounds(:,1)));

%% initialize output
retcode=nan(1,monte_carlo_nsim);
logpost=nan(1,monte_carlo_nsim);
LogLik=nan(1,monte_carlo_nsim);
log_prior=nan(2,monte_carlo_nsim);
behave=false(1,monte_carlo_nsim);
could_solve=false(1,monte_carlo_nsim);
%% find the behaving and solving parameters
for isim=1:monte_carlo_nsim
    obj=assign_estimates(obj,xparam(:,isim));
	[obj,retcode(isim)]=solve(obj);
    could_solve(isim)=~retcode(isim);
    if could_solve(isim) && smooth_flag
        [logpost(isim),LogLik(isim),log_prior(:,isim),junk,retcode(isim),obj]=log_posterior_kernel(obj,xparam(:,isim),smooth_flag);
		
        % I could recover the one-step-ahead predictions here to map the
        % fit of the model... alternatively, this could be done in the next
        % step (monte_carlo_behavioral_function)
    end
    
    if behavior
        % there is an issue here as the behavioral function may or may not
        % require filtering...
        if ~retcode(isim)
            behave(isim)=monte_carlo_behavioral_function(obj,varargs{:});
        end
    else
        behave(isim)=could_solve(isim);
    end
end

disp([num2str(100*sum(could_solve)/monte_carlo_nsim),' percent of the prior support has a solution'])
disp([num2str(100*sum(behave)/monte_carlo_nsim),' percent of the prior support is consistent with the behavior'])

