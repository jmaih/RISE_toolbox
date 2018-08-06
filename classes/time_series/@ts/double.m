function d=double(self)
% DOUBLE: Returns the underlying data of the time series
%
% ::
%
%    data = double(db);
%
% Args:
%
%    db (ts object): time series object
%
% Returns:
%    :
%    - data (numeric): vector/matrix/tensor form of the data underlying the time series
%

d=self.data;
end