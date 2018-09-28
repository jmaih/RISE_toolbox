function out=irf(obj,varargin)
% Computes impulse responses for a RISE model
%
% ::
%
%   myirfs=irf(obj)
%
%   myirfs=irf(obj,varargin)
%
% Args:
%
%    obj (rise | dsge): single or vector of RISE models
%
%    varargin : optional options coming in pairs. The notable ones that
%       will influence the behavior of the impulse responses are:
%
%       - **irf_shock_list** [char | cellstr | {''}]: list of shocks for which we
%         want to compute impulse responses
%
%       - **irf_var_list** [char | cellstr | {''}]: list of the endogenous variables
%         we want to report
%
%       - **irf_periods** [integer | {40}]: length of the irfs
%
%       - **irf_shock_sign** [numeric | -1 | {1}]: sign or scale of the original
%         impulse. If **irf_shock_sign** >0, we get impulse responses to a
%         positive shock. If **irf_shock_sign** <0, the responses are negative.
%         If If **irf_shock_sign** =0, all the responses are 0.
%
%       - **irf_draws** [integer | {50}]: number of draws used in the simulation
%         impulse responses in a nonlinear model. A nonlinear model is defined as
%         a model that satisfies at least one of the following criteria
%
%         - solved at an order >1
%         - has more than one regime and option **irf_regime_specific** below is
%           set to false
%
%       - **irf_type** [{irf} | girf]: type of irfs. If the type is irf, the
%         impulse responses are computed directly exploiting the fact that the
%         model is linear. If the type is girf, the formula for the generalized
%         impulse responses is used: the irf is defined as the expectation of the
%         difference of two simulation paths. In the first path the initial
%         impulse for the shock of interest is nonzero while it is zero for the
%         second path. All other shocks are the same for both paths in a given
%         simulation.
%
%       - **irf_regime_specific** [{true} | false]: In a switching model, we may or
%         may not want to compute impulse responses specific to each regime.
%
%       - **irf_use_historical_data** [{false} | true]: if true, the data stored in
%         option **simul_historical_data** are used as initial conditions. But
%         the model has to be nonlinear otherwise the initial conditions are set
%         to zero. This option gives the flexibility to set the initial
%         conditions for the impulse responses.
%
%       - **irf_to_time_series** [{true} | false]: If true, the output is in the
%         form of time series. Else it is in the form of a cell containing the
%         information needed to reconstruct the time series.
%
% Returns:
%    :
%
%    - **myirfs** [{struct}|cell]: Impulse response data
%
% Note:
%
%    - for linear models or models solved up to first order, the initial
%      conditions as well as the steady states are set to 0 in the computation
%      of the impulse responses.
%
%    - for nonlinear models, the initial conditions is the ergodic mean
%
%    - irf automatically enables simul_bgp_deviation. The intuition is clear:
%      if we apply the GIRFs, the growth components will always fall out.
%      Somehow it is a different story with Tz_sig
%

if isempty(obj)
    
    mydefaults=irf@generic(obj);
    
    mydefaults=[mydefaults
        {'irf_anticipate',true,@(x)islogical(x),...
        'irf_anticipate must be true or false'}];
        
    if nargout
        
        out=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end

    
    return
    
end

obj=set(obj,varargin{:});

for iobj=1:numel(obj)
    
    obj(iobj).options.simul_anticipate_zero=~obj(iobj).options.irf_anticipate;
    
end

out=irf@generic(obj);

end