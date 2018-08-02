function oo_=moments(db,drange,ar,lambda)
% Computes the empirical moments of a time series
%
% ::
%
%    oo_ = moments(db);
%    oo_ = moments(db,drange);
%    oo_ = moments(db,drange,ar);
%    oo_ = moments(db,drange,ar,lambda);
%
% Args:
%
%    db (ts object): time series object to get data
%    drange (char | serial date | cellstr | {[]}): Range of the data to use
%    ar (integer | {1}): order of autocorrelation
%    lambda (numeric | {[]}): hyperparameter for hp-filtering the data
%      before computing the moments. If empty, the data are not
%      hp-filtered.
%
%
% Returns:
%
%    - **oo_** [struct]: structure with fields
%
%       - vcov : variance covariance
%       - skewness : skewness
%       - kurtosis : kurtosis
%       - variance : variance
%       - stdev : standard deviation
%       - corrcoef : correlation array
%

if nargin<4

    lambda=[];

    if nargin<3

        ar=[];

        if nargin<2

            drange=[];

        end

    end

end

if isempty(ar)

    ar=1;

end

if ~isempty(drange)

    db=db(drange);

end

oo_.mean=mean(db);

if ~isempty(lambda)

    [oo_.hp_trend,oo_.hp_cycle] = hpfilter(db,lambda);

    db=oo_.hp_cycle;

else

    db = bsxfun(db,@minus,oo_.mean);

end

oo_.vcov=cov(db);

oo_.skewness=skewness(db);

oo_.kurtosis=kurtosis(db);

oo_.variance = diag(oo_.vcov);oo_.variance=oo_.variance(:).';

oo_.stdev = sqrt(oo_.variance);

oo_.corrcoef = corrcoef(db);

if ar > 0

    y=double(db);

    for i=1:ar

        oo_.autocorr{i} = y(ar+1:end,:)'*y(ar+1-i:end-i,:)./((size(y,1)-ar)*std(y(ar+1:end,:))'*std(y(ar+1-i:end-i,:)));

    end

end

end