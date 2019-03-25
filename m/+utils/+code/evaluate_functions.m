%  INTERNAL FUNCTION: Evaluates functions in various formats
% 
%  ::
% 
%    varargout=evaluate_functions(xcell,varargin)
% 
%  Args:
% 
%     - **xcell** [cell|struct|fhandle]: function of interest
% 
%       - input is a cell array: each element in a cell is assumed to be a
%         function handle
%       - input is a struct:
% 
%         - if the field names are among: 'size','functions','map','partitions'
%           then we are dealing with derivatives
%         - if the field names are among: 'code','argins','argouts' then the
%           function could not be written as a function handle and has to be
%           evaluated using eval
% 
%       - input is a function handle : direct evaluation
% 
%     - **varargin** : input arguments to the function
% 
%  Returns:
%     :
% 
%     - **varargout** : ouput arguments to the function
% 
%