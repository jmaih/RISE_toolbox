function start=update_posterior_simulation_initial_conditions(obj,start,new_vcov,acceptance_rate)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 


nobj=numel(obj);
if nobj==0
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    start=struct();
    return
end

if obj.markov_chains.regimes_number==1 && obj.options.vp_analytical_post_mode
    % do nothing ...
else
    start=update_posterior_simulation_initial_conditions@rise_generic(obj,...
        start,new_vcov,acceptance_rate);
end
