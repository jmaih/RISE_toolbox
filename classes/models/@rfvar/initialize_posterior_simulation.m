function [init,sampler,total_draws]=initialize_posterior_simulation(obj)
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
    init=struct();
    return
end

if obj.markov_chains.regimes_number==1 && obj.options.vp_analytical_post_mode
    % [smpl,fsmpl,accept_rate,start] = constant_bvar_sampler(start,nsamples,varargin)
    sampler=@utils.mcmc.constant_bvar_sampler;
    number_of_burns=round(obj.options.mcmc_burn_rate*obj.options.mcmc_number_of_simulations);
    if isempty(obj.constant_var_data)
        error('The posterior mode needs to be found first')
    end
    init=obj.constant_var_data;
    total_draws=number_of_burns+obj.options.mcmc_number_of_simulations;
else
    [init,sampler,total_draws]=initialize_posterior_simulation@rise_generic(obj);
end

