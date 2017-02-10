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
% - Keeps the start date and only changes the data. New data no longer
% needs to be of the same size as the old one.
%
% Examples
% ---------
%
% See also:

ts.check_size(newdata);

this.data=newdata;

end
