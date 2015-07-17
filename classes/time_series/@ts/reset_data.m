function this=reset_data(this,newdata)
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
if ~isequal(size(this.data),size(newdata))
    error('data size mismatch')
end
this.data=newdata;
end
