function varargout=resort(varargin)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

nout=nargout;
if nout~=nargin
    error([mfilename,':: number of inputs and outputs should be the same'])
end
varargout=cell(1,nout);
[varargout{1},tags]=sort(varargin{1});
for ii=2:nout
    varargout{ii}=varargin{ii}(:,tags);
end

end
