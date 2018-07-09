function varargout=decipher(varargin)
% INTERNAL FUNCTION
%

% Interprets error codes
%
% ::
%
%   msgout=decipher(code)
%
% Args:
%
%    code (scalar | vector): scalar or vector of error codes returned by RISE.
%
% Returns:
%    :
%
%    - **msgout** [char|cellstr]: explanation of the return code as a char if
%      the input is a scalar or as a cellstr if input is a vector.
%

[varargout{1:nargout}]=utils.error.decipher(varargin{:});

end
