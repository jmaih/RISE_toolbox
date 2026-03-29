%--- help for generic/theoretical_autocorrelations ---
%
%  theoretical_autocorrelations computes (mostly) regime-specific
%  partial autocorrelations
% 
%  Args:
% 
%    - obj : [rise|dsge] : model object
%    - varargin : pairwise arguments with the most relevant for this
%      function being : 
% 
%      - autocorr_ar : [numeric>=0|{5}] : order of autocorrelation. If you want to
%        compute the contemporaneous correlation only set to order to 0
%      - autocov_aggregate : [true|{false}] : when true and under regime
%        switching, the autocovariances are aggregated across regimes with the
%        weights being the ergodic distribution of the regimes. Then the
%        autocorrelations are computed out of the aggregated autocovariances.
%        Otherwise, by default, the theoretical autocorrelations are regime-specific
%      - autocov_model_resolve : [{true}|false] : if true, the function tries
%        to solve/resolve the model first    
% 
%  Returns:
%     :
% 
%     - **Acorr** [cell array]: cell array of regime-specific correlations, 
%       when autocov_aggregate is false. When there is only one regime or when
%       autocov_aggregate is true, the cell array contains only one element.
%       Inside each element is Acorr is a possibly 3-dimensional array where the
%       3rd dimension represents the order of autocorrelation. In order words,
%       Acorr{regime=i}(:,:,1) is the contemporaneous correlation while
%       Acorr{regime=i}(:,:,2) is the autocorrelation of order 1, so that
%       Acorr{regime=i}(:,:,k+1) is the autocorrelation of order k.
% 
%     - **Acov** [cell array]: cell array of regime-specific autocovariances, 
%       following the same logic as the autocorrelations described above
% 
%     - **retcode** [scalar]: return code
% 
%  Note:
% 
%     - The sucessful computation of the autocorrelations depends on the
%       solving of the discrete Lyapunov or Stein equation.
%