function a=anonymize_to_one_input(main,varargin)
% anonymize_to_one_input - creates an anonymous function with one input
%
% Syntax
% -------
% ::
%
%   a=anonymize_to_one_input(main,varargin)
%
% Inputs
% -------
%
% - **main** [function handle]: function handle taking n+1 arguments
%
% - **varargin** [undefined]: arguments #2 to #n+1 of function "main"
%
% Outputs
% --------
%
% - **a** [function handle]: anonymous function
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

a=@(y)main(y,varargin{:});
end