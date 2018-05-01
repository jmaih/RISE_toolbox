function this=unfold(last_hist_date,varargin)
% UNFOLD -- takes a folded times series and turns it into a one-page time
% series
%
% ::
%
%
%   this=UNFOLD(last_hist_date,obj)
%
%   this=UNFOLD(last_hist_date,obj1,obj2,...,objn)
%
% Args:
%
%    - **last_hist_date** [char|serial date]:
%
%    - **obj** [ts]: time series
%
% Returns:
%    :
%
%    - **this** [ts]: unfolded time series
%
% Note:
%
% Example:
%
%    See also:

if ~(ischar(last_hist_date)||isnumeric(last_hist_date))
    error('pivot date should be a valid date or a serial date')
end

obj=ts.collect(varargin{:});
sdate=date2serial(last_hist_date);
dn=obj.date_numbers;
loc=find(dn==sdate);
if isempty(loc)
    error('pivot date not found')
end
data=double(obj);
first_page=data(1:loc,:,1);
extras=squeeze(data(loc,:,2:end));
first_page=[first_page;extras.'];
this=ts(obj.start,first_page,obj.varnames);

end