function hdl=psrf_plot(obj,pname)
% Make a plot of the "posterior scale reduction factor," i.e., Gelman-Rubin diagnotics from the chains
%
% ::
%
%    hdl = psrf_plot(obj, pname)
%
% Args:
%    obj (mcmc object): mcmc object
%    pname (str): parameter name
%
% Returns:
%    :
%
%    - **hdl** (handle object): handle to plot object
%
% Warning:
%    - This function requires multiple chains of MCMC samples. See
%      **nchain** option of samplers. 
%
% References:
%    - :cite:`gelman1992inference`
%

% Reference:
%    - Gelman, Andrew, and Donald B. Rubin. "Inference from iterative
%      simulation using multiple sequences." Statistical science 7.4 (1992):
%      457-472.  
%
hdl0=plot(obj.psrf.(pname).time,obj.psrf.(pname).psrf);

title(sprintf('%s(%0.4f)',pname,obj.psrf.(pname).psrf(end)))

axis tight

if nargout

    hdl=hdl0;

end

end