function fname=remove_file_extension(fname)
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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

thedot=find(fname=='.');
if ~isempty(thedot)
    fname=fname(1:thedot-1);
end
end