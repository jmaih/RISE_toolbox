function [obj,H,issue]=hessian(obj,varargin)
% Computes the hessian of the model at a specific point
%
% ::
%
%   [obj,H,issue] = hessian(obj)
%
%   [obj,H,issue] = hessian(obj,x)
%
%   [obj,H,issue] = hessian(obj,x,varargin)
%
% Args:
%
%    obj (rise | dsge | rfvar | svar): model object
%
%    x ([] | vector): vector at which one wants to compute the hessian
%
%    varargin: additional optional inputs among which the most relevant
%      for estimation is:
%
%      - **hessian_type** [{'fd'}|'opg']: The hessian is either computed by
%        finite differences (fd) or by outer-product-gradient (opg)
%
% Returns:
%    :
%
%    - **obj** [rise|dsge|rfvar|svar]: model object containing the new hessian
%      in case the model was previously estimated.
%
%    - **H** [d x d matrix]: hessian matrix
%
%    - **issue** [char|'']: any issue encountered during the computation of
%      the hessian
%

nobj=numel(obj);

if nobj==0

    if nargout>1

        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])

    end

    obj=utils.hessian.numerical();

    return
    
elseif nobj>1
        
    NumWorkers=utils.parallel.get_number_of_workers();
        
    if NumWorkers
        
        pctRunOnAll('utils.estim.warnings_disable();')
        
        pctRunOnAll('utils.estim.prior.warnings_fmincon_disable();')
        
    end
    
    parfor (ii=1:nobj,NumWorkers)
        
        obj(ii)=hessian(obj(ii),varargin{:}); %#ok<PFBNS>
        
    end
    
    return

end

x=[];

if nargin>1
    
    x=varargin{1};
    
    varargin=varargin(2:end);
    
end    

if isempty(x)

    if ~isfield(obj.estimation.posterior_maximization,'mode')||...
            isempty(obj.estimation.posterior_maximization.mode)

        error('vector of parameters should be specified when model has not been estimated')

    end

   x=obj.estimation.posterior_maximization.mode;

end

if ~isempty(varargin)

    obj=set(obj,varargin{:});

end

ws=utils.estim.warnings_disable();

fh=pull_objective(obj);

[H,issue]=utils.hessian.numerical(fh,x,obj.options.hessian_type);

utils.estim.warnings_enable(ws);

post_max=obj.estimation.posterior_maximization;

post_max.hessian=H;

post_max=generic_tools.posterior_maximization_variable_quantities(post_max,...
    obj.linear_restrictions_data.a_func);

obj.estimation.posterior_maximization=post_max;

end