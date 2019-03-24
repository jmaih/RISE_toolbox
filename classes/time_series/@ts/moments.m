%--- help for ts/moments ---
%
%  Computes the empirical moments of a time series
% 
%  ::
% 
%     oo_ = moments(db);
%     oo_ = moments(db,drange);
%     oo_ = moments(db,drange,ar);
%     oo_ = moments(db,drange,ar,lambda);
% 
%  Args:
% 
%     db (ts object): time series object to get data
% 
%     drange (char | serial date | cellstr | {[]}): Range of the data to use
% 
%     ar (integer | {1}): order of autocorrelation
% 
%     lambda (numeric | {[]}): hyperparameter for hp-filtering the data
%       before computing the moments. If empty, the data are not
%       hp-filtered.
% 
% 
%  Returns:
% 
%     - **oo_** [struct]: structure with fields
% 
%        - vcov : variance covariance
%        - skewness : skewness
%        - kurtosis : kurtosis
%        - variance : variance
%        - stdev : standard deviation
%        - corrcoef : correlation array
% 
%
%    Other functions named moments
%
%       mdd/moments
%