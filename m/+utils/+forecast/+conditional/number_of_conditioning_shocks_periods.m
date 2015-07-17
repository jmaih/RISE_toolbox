function ncsp = number_of_conditioning_shocks_periods(hypo,ncp,nap)
% number_of_conditioning_shocks_periods - length of shocks over the forecast horizon
%
% Syntax
% -------
% ::
%
%   ncsp = number_of_conditioning_shocks_periods(hypo,ncp,nap)
%
% Inputs
% -------
%
% - **ncp** [integer]: number of periods over which we have conditioning
%   information
%
% - **nap** [integer|{size(G,3)}]: number of anticipated periods i.e. how
%  far agents see into the future + the current period
%
% - **hypo** [NCP|NAS|JMA]: forecasting hypothesis, determining the
%   number of periods over which future shocks will be drawn or assumed
%   known to the agents.
%   - **NCP** : under this scheme, the number of periods of future shocks 
%       is equal to the number of periods over which conditioning
%       information is available. Beyond the conditioning period, the
%       shocks are 0.
%   - **NAS** : under this scheme, the number of periods of future shocks 
%       is equal to the number of anticipated steps. Beyond that number,
%       the shocks are expected to be 0.
%   - **JMA** : under this scheme, the number of periods of future shocks 
%       is equal to the number of anticipated steps plus the number of
%       conditioning periods. It is assumed that every condition is
%       determined by the same number of future shocks
%
% Outputs
% --------
%
% - **ncsp** [integer]: number of periods over which future shocks are
%   computed 
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

hypo=upper(hypo);
switch hypo
    case {'NCP'}
        ncsp=ncp;
    case {'NAS'}
        ncsp=nap;
        if ~isequal(ncp,nap)
            error([mfilename,':: for the NAS assumption, you need # anticipated steps = # conditioning periods'])
        end
    case {0,'JMA'}
        ncsp=nap+ncp;
    otherwise
        error([mfilename,':: Unknown option for the anticipation hypothesis'])
end

end
