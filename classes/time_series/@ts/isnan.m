function flag=isnan(self)
% Returns whether the corresponding data is nan
%
% ::
%
%    flag = isnan(db);
%
% Args:
%
%    db (ts object): time series object
%
% Returns:
%    :
%
%    - flag (bool): whether the corresponding data is nan or not
%
% Note:
%    - Since the time series object supports logical indexing, one can
%      directly use the resulting flag with the time series object. 
%

flag=isnan(self.data);
end