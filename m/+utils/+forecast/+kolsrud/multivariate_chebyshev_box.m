function [mvcb,my]=multivariate_chebyshev_box(y,gam,c)
% multivariate_chebyshev_box - constructs chebyshev boxes for
% multivariate-multiperiods densities
%
% ::
%
%
%   mvcb=multivariate_chebyshev_box(y,gam,c)
%
% Args:
%
%    - **y** [numeric] : N x T x G array, with
%      - **N** [numeric] : number of simulations/replications
%      - **T** [numeric] : sample length (time series dimension)
%      - **G** [numeric] : number of variables
%
%    - **gam** [scalar|vector] : percentile(s)
%
%    - **c** [empty|scalar|vector] : precomputed chebyshev distances
%
% Returns:
%    :
%
%    - **mvcb** [2 x T x G x numel(gam)] :  array of boxes
%
%    - **my** [1 x T x G] :  mean across simulations
%
% Note:
%
% Example:
%
%    See also: chebyshev_distance, standardized_distance

% References:
% Dag Kolsrud (2015): "A time-simultaneous prediction box for a
% multivariate time series", Journal of Forecasting

if nargin<3
    c=utils.forecast.kolsrud.chebyshev_distance(y);
end
N=numel(c);
[~,tag]=sort(c);
y=y(tag,:,:);
my=mean(y,1);

if ~issorted(gam)
    error('gam must be sorted')
end
mvcb=[];
for ig=1:numel(gam)
    do_one_box()
end

    function do_one_box()
        cutoff=ceil(N*gam(ig));
        yM=y(1:cutoff,:,:);
        if ig==1
            mvcb=make_mvcb_data();
        else
            mvcb(:,:,:,ig)=make_mvcb_data();
        end
        function d=make_mvcb_data()
            d=[
                min(yM,[],1)
                max(yM,[],1)
                ];
        end
    end
end