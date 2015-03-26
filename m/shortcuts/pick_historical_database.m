function varargout=pick_historical_database(varargin)
% pick_historical_database pick a historical database in the filters
%
% Syntax
% -------
% ::
%
%   histdb=pick_historical_database(m)
%
%   histdb=pick_historical_database(m,type)
%
%   histdb=pick_historical_database(m,type,hist_end_date)
%
% Inputs
% -------
%
% - **m** [rise|dsge] : model object
%
% - **type** [{'mean'}|'random'] : type of database. If 'mean', the
% database with expected values is returned. If 'random', a database is
% drawn randomly using the probabilitiy distribution at the end of history.
%
% - **hist_end_date** [char|{sample end date}] : desired date for the end
% of history
%
% Outputs
% --------
%
% - **histdb** [struct] : structure with the time series of the endogenous
% variables, the regime probabilities
%
% More About
% ------------
%
% - if the number of regimes is 1, "type" and "hist_end_date" are
% superfluous.
%
% - if a database is picked at random, the regime probabilities are changed
% so as to inform the forecasting procedure later on about what regime the
% data belong to. For instance, suppose we have 3 regimes and that the
% second regime is picked based on the distribution [.5,.2,.3] at the
% desired (or default) end date of history. Once the regime is picked,
% those probabilities are modified to [0,1,0].
%
% Examples
% ---------
%
% See also: 

[varargout{1:nargout}]=utils.miscellaneous.pick_historical_database(varargin{:});

end