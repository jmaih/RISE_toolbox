%--- help for mcmc/plot_brooks_gelman ---
%
%  Make a plot of cumulative means from the mcmc chain
% 
%  ::
% 
%     h=plot_brooks_gelman(this,myoutput)
%     h=plot_brooks_gelman(this,myoutput,nparamsPerFigure)
% 
%  Args:
% 
%     this (mcmc object): mcmc object
% 
%     myoutput (struct): output of brooks_gelman
% 
%     nparamsPerFigure [numeric|{3}]: number of univariate parameters per
%     figure
% 
%  Returns:
%     :
% 
%     - **h** (handle object): handle to all the plots (unvariates and multivariate)
%