%  pick_historical_database pick a historical database in the filters
% 
%  ::
% 
% 
%    histdb=pick_historical_database(m)
% 
%    histdb=pick_historical_database(m,type)
% 
%    histdb=pick_historical_database(m,type,hist_end_date)
% 
%  Args:
% 
%     - **m** [rise|dsge] : model object
% 
%     - **type** [{'mean'}|'random'] : type of database. If 'mean', the
%     database with expected values is returned. If 'random', a database is
%     drawn randomly using the probabilitiy distribution at the end of history.
% 
%     - **hist_end_date** [char|{sample end date}] : desired date for the end
%     of history
% 
%  Returns:
%     :
% 
%     - **histdb** [struct] : structure with the time series of the endogenous
%     variables, the regime probabilities
% 
%  Note:
% 
%     - The database is always drawn from the updated data and not the smoothed
%     ones.
% 
%     - if the number of regimes is 1, "type" and "hist_end_date" are
%     superfluous.
% 
%     - if a database is picked at random, the regime probabilities are changed
%     so as to inform the forecasting procedure later on about what regime the
%     data belong to. For instance, suppose we have 3 regimes and that the
%     second regime is picked based on the distribution [.5,.2,.3] at the
%     desired (or default) end date of history. Once the regime is picked,
%     those probabilities are modified to [0,1,0]. This information is useful
%     for determining the distribution, the next regime is drawn from.
% 
%     - if the type of database chosen is the "mean", the probability
%     distribution of the historical regimes remain unchanged. i.e. [.5,.2,.3]
%     remains [.5,.2,.3]. The consequence is that the historical state is not
%     known with certainty.
% 
%     - if the user adds a time series called "regime" to the database, the
%     information on the distribution of the historical probabilities is
%     irrelevant since the regime variables should inform the forecasting
%     procedure about where the system is going in the future.
% 
%  Example:
% 
%     See also:
%