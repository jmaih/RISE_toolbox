function flag=isinf(self)
% Returns whether the corresponding data is infinite
%
% ::
%
%    flag = isinfinite(db);
%
% Args:
%    db (ts object): time series object
%
% Returns:
%    :
%
%    - flag (bool): whether the corresponding data is infinite or not
%
% Note:
%    - Since the time series object supports logical indexing, one can directly use the resulting flag with the time series object.
%

flag=isinf(self.data);
end