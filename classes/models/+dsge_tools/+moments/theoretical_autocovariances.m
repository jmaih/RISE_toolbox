%--- help for generic/theoretical_autocovariances ---
%
%  `THEORETICAL_AUTOCOVARIANCES` - Compute regime-specific autocovariances
% 
%  Syntax:
%    [V,retcode] = theoretical_autocovariances(obj)
%    [V,retcode] = theoretical_autocovariances(obj, autocov_ar)
%    [V,retcode] = theoretical_autocovariances(obj, autocov_ar, vList)
%    [V,retcode] = theoretical_autocovariances(obj, autocov_ar, vList, varargin)
% 
%  Args:
%    - `obj` : [rise|dsge] : Model object
%    - `autocov_ar` : [numeric>=0|{5}] : Order of autocovariance. If you want to
%      compute the variance only, set the order to 0.
%    - `vList` : [empty|cellstr|vector] : Names of indices/positions for the
%      variables for which to compute the autocovariance.
%    - `varargin` : Additional options specified as key-value pairs.
% 
%  Returns:
%    - `V` : [cell array] : Cell array of regime-specific covariances. When
%      autocov_aggregate is false, each element of V is a 3-dimensional array,
%      where the 3rd dimension represents the order of autocovariance. In other
%      words, V{regime=i}(:,:,1) is the variance, and V{regime=i}(:,:,k+1) is the
%      autocovariance of order k.
%    - `retcode` : [scalar] : Return code. 0 indicates success.
% 
%  Note:
%    - The successful computation of the autocovariances depends on the
%      solution of the discrete Lyapunov or Stein equation. Different Lyapunov
%      algorithms can be set through the option "lyapunov_algo". The available
%      algorithms are:
% 
%        - 'doubling'  OR @lyap_solvers.doubling (the default)
%        - 'bicgstab'  OR @lyap_solvers.bicgstab
%        - 'bicgstabl' OR @lyap_solvers.bicgstabl
%        - 'cgs'       OR @lyap_solvers.cgs
%        - 'fevdcov'   OR @lyap_solvers.fevdcov
%        - 'direct'    OR @lyap_solvers.direct
%        - 'fix_point' OR @lyap_solvers.fix_point
%        - 'robust'    OR @lyap_solvers.robust
%        - 'sandwich'  OR @lyap_solvers.sandwich
%        - 'schur'     OR @lyap_solvers.schur
%        - 'bartels_stewart' OR @lyap_solvers.bartels_stewart
% 
%    - In the presence of a nonstationary model, it is not a good idea to use
%      the doubling algorithm as it won't converge. A different algorithm, such
%      as "bartels_stewart," may work better.
% 
%    - Options autocov_aggregate, autocov_aggregate, and autocov_ar can be
%      passed as pairs in varargin. In the case of autocov_ar, passing it
%      directly overrides the version passed in varargin.
% 
%  See also:
%    - `SOLVE`, `LYAPUNOV_EQUATION`
%