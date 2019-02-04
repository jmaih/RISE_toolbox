function hdl=psrf_plot(obj,pname,start)
% Make a plot of the "posterior scale reduction factor," i.e., Gelman-Rubin
% diagnotics from the chains 
%
% ::
%
%    hdl = psrf_plot(obj, pname)
%
% Args:
%
%    obj (mcmc object): mcmc object
%
%    pname (char): parameter name. N.B: One of the parameter names is
%    "multivariate_" and it represents the aggregated statistics.
%
%    start (numeric|{1}|function handle): iteration at which to start the
%       plot of the PSRF. If a function handle is used then it should take
%       as input the total number of observations and return the point at
%       which to start. e.g. @(x)round(0.5*x)
%
%
% Returns:
%    :
%
%    - **hdl** (handle object): handle to plot object
%
% Warning:
%
%    - This function requires multiple chains of MCMC samples. See
%      **nchain** option of samplers. 
%
% References:
%
%    - :cite:`gelman1992inference`
% 

% Reference:
%    - Gelman, Andrew, and Donald B. Rubin. "Inference from iterative
%      simulation using multiple sequences." Statistical science 7.4 (1992):
%      457-472.  
%

if nargin<3
    
    start=[];
    
end

if isempty(start)
    
    start=1;
    
elseif isnumeric(start)
    
    if ~isscalar(start)||...
            ~isfinite(start)||...
            ~(floor(start)==ceil(start))||...
            start<1
        
        error('start must be a scalar finite and positive integeger')
        
    end
    
elseif isa(start,'function_handle')
    
    n=numel(obj.psrf.(pname).time);
    
    start=start(n);
        
    
end

hdl0=plot(obj.psrf.(pname).time(start:end),obj.psrf.(pname).psrf(start:end));

title(sprintf('%s(%0.4f)',pname,obj.psrf.(pname).psrf(end)))

axis tight

if nargout

    hdl=hdl0;

end

end