function this=fold(last_hist_date,varargin)
% FOLD -- takes the first page of possibly several databases and adds
% further pages using the information from a pivot date
%
% ::
%
%
%   this=FOLD(last_hist_date,obj)
%
%   this=FOLD(last_hist_date,obj1,obj2,...,objn)
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
%    - **this** [ts]: folded time series
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
nrows=size(first_page,1);
extras=data(loc+1:end,:);
[nx,ncols]=size(extras);
extra_pages=nan(nrows,ncols,nx);

extra_pages(end,:,:)=extras.';
first_page=cat(3,first_page,extra_pages);

this=ts(obj.start,first_page,obj.varnames);
end