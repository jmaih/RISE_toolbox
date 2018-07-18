function [C,R]=autocov(self,params,max_periods)
% Compute the autocovariances (and the auto-correlation) of endogenous variables given the parameter values
%
% ::
%
%    [C,R] = autocov(self);
%    [C,R] = autocov(self, params);
%    [C,R] = autocov(self, params, max_periods);
%
% Args:
%    self (var object): var object
%    params (cell of struct): struct containing var model related parameters (default: [])
%    max_periods (integer): maximum number of period to calculate auto-covariance (default: 5)
%
% Returns:
%    :
%
%    - **C** [4-dimensional array]: auto-covariance where dimensions correspond to
%
%       - 1,2: covariance
%       - 3: time lags
%       - 4: Number of parameters
%
%    - **R** [4-dimensional array]: auto-correlation with same dimensions
%

if nargin<3

    max_periods=[];

    if nargin<2

        params=[];

    end

end

if isempty(max_periods),max_periods=5; end

params=solve(self,params);

Rfunc=identification(self,'choleski');

np=numel(params);

for ip=1:np

    if ip==1

        [C,R,info]=do_one_parameter(params(:,ip)); %#ok<ASGLU>

        C=C(:,:,:,ones(1,np));

        R=R(:,:,:,ones(1,np));

    else

        [C(:,:,:,ip),R(:,:,:,ip)]=do_one_parameter(params(:,ip));

    end

end

    function [C,R,info]=do_one_parameter(param)

        [C,R,info]=vartools.autocorr(param.B,...
            Rfunc(param),...
            self.nx*self.ng,max_periods);

    end

end
