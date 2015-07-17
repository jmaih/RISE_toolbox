function [this1,this2]=intersect(this1,this2)
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

if nargin~=2
    error([mfilename,':: function must have 2 arguments'])
end
if ~isa(this1,'ts')||~isa(this2,'ts')
    error([mfilename,':: both arguments should be time series'])
end
if ~strcmp(this1.frequency,this2.frequency)
    error([mfilename,':: datasets must have same frequency'])
end
if this1.NumberOfPages~=this2.NumberOfPages
    error([mfilename,':: datasets must have same number of pages'])
end
[~,I1,I2] = intersect(this1.date_numbers,this2.date_numbers);
if isempty(I1)||isempty(I2)
    error([mfilename,':: arguments must have common dates'])
end
C=this1.date_numbers(I1);
vname1=this1.varnames;
vname2=this2.varnames;
this1=ts(C,this1.data(I1,:,:),vname1);
this2=ts(C,this2.data(I2,:,:),vname2);
end
