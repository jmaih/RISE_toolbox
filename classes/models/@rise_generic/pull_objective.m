function [ff,lb,ub]=pull_objective(obj,varargin)
% pull_objective -- pulls the objective function to optimize
%
% Syntax
% -------
% ::
%
%   [ff,lb,ub]=pull_objective(obj)
%
%   [ff,lb,ub]=pull_objective(obj,varargin)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|svar|rfvar]: initial model object
%
% - **varargin** [pairwise addional inputs]: usual RISE arguments
%
% Outputs
% --------
%
% - **ff** [function handle]: for computing "minus log posterior kernel"
%
% - **lb** [d x 1 vector]: lower bound of the parameters to optimize
%
% - **ub** [d x 1 vector]: upper bound of the parameters to optimize
%
% More About
% ------------
%
% - The function can be used for :
%   - optimization, 
%   - gradient computation, 
%   - hessian computation, 
%   - posterior simulation 
%
% - Using this function is potentially costly, one could alternatively
% simply use log_posterior_kernel. However, if there are restrictions, they
% will not be enforced. Nevertheless it is an interesting proposition that
% should be investigated further.
%
% Examples
% ---------
%
% See also: 

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    ff=struct();
    return
end
lb=[obj.estimation.priors.lower_bound]';
ub=[obj.estimation.priors.upper_bound]';

% try to avoid unnecessary computations like storing filters and so on
%----------------------------------------------------------------------
obj.estimation_under_way=true;
obj.options.kf_filtering_level=0;

ff=@engine;

    function [minus_log_post,retcode]=engine(x)
        retcode=0;
        if any(x<lb)||any(x>ub)
            retcode=7;
        else
            [~,minus_log_post,~,issue,viol]=estimation_wrapper(obj,[],x,lb,ub,0);
            % function evaluations computed elsewhere
            if ~isempty(issue)
                retcode=issue;
            end
            if ~isempty(viol) && any(viol>0)
                retcode=7;
            end
        end
        if retcode
            minus_log_post=obj.options.estim_penalty;
        end
    end
end