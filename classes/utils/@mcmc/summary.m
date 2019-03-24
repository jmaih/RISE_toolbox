%--- help for mcmc/summary ---
%
%  Print summary information about the mcmc draws
% 
%  ::
% 
%     [myMeanStdev, MyQuantiles] = summary(mcobj);
%     [myMeanStdev, MyQuantiles] = summary(mcobj, varargin);
% 
%  Args:
% 
%     mcobj (mcmc object): mcmc object
% 
%     varargin (options): options need to come in pairs:
% 
%        - 'percnt': Quantile points (default: [2.5, 25, 50, 75, 97.5])
%        - 'batch_size': sampling size of for the computation of thinned standard deviation (SD(batch)), every 'batch_size'th points are used to compute SD(batch).
% 
%  Returns:
%     :
% 
%     - **myMeanStdev** [cell]: mean and standard deviations
%     - **myQuantiles** [cell]: quantile values
% 
%
%    Other functions named summary
%
%       categorical/summary      dataset/summary        tabular/summary
%       codistributed/summary    distributed/summary    tall/summary
%