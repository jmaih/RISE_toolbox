%--- help for abstvar/setOptions ---
%
%  INTERNAL FUNCTION: Sets options for the VAR
% 
%  ::
% 
%     var_obj = set_inputs(var_obj,varargin);
% 
%  Args:
%     var_obj (var object): VAR object
%     varargin: Options to set. Must come in pairs
% 
%        - **'data'**: Data for the VAR. Can be either a struct or ts object
%        - **'prior'**: Prior to use if BVAR. Available priors are
% 
%           - minnesota
%           - jeffrey
%           - nw
%           - normal-wishart
%           - inw
%           - sz
% 
%        - **'linear_restrictions'**: Linear restrictions. For the format, see XXXXXXXXX
% 
%  Returns:
%     :
% 
%        - var_obj (var object): New VAR with the given inputs updated.
% 
%
%    Other functions named setOptions
%
%       AbstractPortfolio/setOptions    PortfolioCVaR/setOptions
%       Portfolio/setOptions            PortfolioMAD/setOptions
%