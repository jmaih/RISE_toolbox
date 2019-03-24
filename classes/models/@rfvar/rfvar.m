%  Reduced Form VAR object
% 
%  This is a class that helps with reduced form VAR analysis. For more detailed explanations on VAR analysis, refer to the `documentation <https://sehyoun.com>`_ .
% 
%  Attributes:
%     endogenous (cellstring): cell of endogenous variable names
%     exogenous (cellstring): cell of exogenous variable names
%     parameters (cellstring): INTERNAL STATES: cell of parameters names
%     nonvar_parameters (cellstring): Non-VAR parameters (depends on markov switching probability)
%     members :XXXXXXXXXXX
%     constant (bool): whether to include constants
%     homogeneity (string): XXXXXXXXXXXX
%     debug (bool): XXXXXXXXX
%     markov_chains ():XXXXXXXXXX
%     is_time_varying_trans_prob ():XXXXXXXXX
%     is_switching (bool): whether a markov switching VAR
%     is_panel (bool): whether a penal VAR or not
%     nlags (integer): number of lags to include
%     nx (integer): number of exogenous variables (constants count as exogenous variables)
%     ng (integer): number of "countries," i.e., the dimensionality of the penal
%     nvars (integer): number of endogenous variables
%     nparams (integer): number of parameters
%     nregs (integer): XXXXX
%     panel_types (cellstring): available panel types
%     panel_with_constant (cellstring): available penal types with constants
% 
%
%    Reference page in Doc Center
%       doc rfvar
%
%