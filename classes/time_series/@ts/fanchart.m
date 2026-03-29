%--- help for ts/fanchart ---
%
%  FANCHART Creates data for plotting a fan chart.
% 
%  out = fanchart(this, ci) generates the data required to plot a fan
%  chart. The input "this" is a time series object with either multiple
%  pages or multiple columns, but not both. The input "ci" is a vector
%  specifying the confidence levels in the range [0, 1] or [0, 100] to be
%  used in calculating the width of the fans.
% 
%  The output "out" is a structure with the following fields:
% 
%  - ci: Vector of confidence levels (same as input ci)
% 
%  - median: Vector of medians of the data
% 
%  - variance: Vector of variances of the data
% 
%  - quantiles: Matrix of quantiles defined by ci
% 
%  - prob_index: Vector of locations of the cutoffs in the data
% 
%  - probs: Vector of probabilities for the quantiles (function of ci)
% 
%  - dates: Vector of serial dates for reconstructing the time series
% 
%  Example:
% 
%  - Generate data for a fan chart and plot it:
% 
%  this = ts('1990Q2', rand(100, 1000));
% 
%  out = fanchart(this, [30, 50, 70, 90]);
% 
%  plot_fanchart(out, 'r', 10);
% 
%  Note:
% 
%  - If the input time series contains several variables, the output is a
%    structure with the names of the different variables in the first level
%    of the fields.
% 
%  See also PLOT_FANCHART.
%