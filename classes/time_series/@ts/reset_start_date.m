function this=reset_start_date(this,startdate)
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


this=ts(startdate,this.data,this.varnames);
end
