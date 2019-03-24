%--- help for ts/spectrum ---
%
%  Computes the spectral density of the data
% 
%  ::
% 
%     [sw,jj,T]=spectrum(this)
% 
%  Args:
% 
%     this (ts object): time series object
% 
%  Returns:
%     :
%     - **sw** (matrix|vector): spectrum of potentially multiple time series
%       in columns
%     - **jj** (vector): range of the spectrum (x-axis)
%     - **T** (integer): number of observations
% 
%