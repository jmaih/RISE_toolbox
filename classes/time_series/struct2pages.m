%  INTERNAL FUNCTION: Brings together several time series objects into a one time series
% 
%  ::
% 
%     this=struct2pages(v1,v2,...,vn)
% 
%  Args:
% 
%     - **Vi** [cell|ts|struct]: time series in ts format:
%     - cell: When Vi is a cell, then its format should be {vname,ts} i.e.
%       the first element is the name of the variable and the second is the
%       data for the variable. In this case, the data must be a single time
%       series
%     - ts:
%     - struct: the fields of the structure should be of the ts format.
% 
%  Returns:
%     :
% 
%     - **this** [ts]: a time series with many columns and potentially many
%       pages
% 
%  See also:
%     ts.collect pages2struct
% 
%