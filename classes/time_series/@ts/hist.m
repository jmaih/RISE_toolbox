function h=hist(varargin)
% Make histograms from time series
%
% ::
%
%    h = hist(db);
%    h = hist(xrange, db);
%    h = hist(xrange, db, M);
%    h = hist(xrange, db, X);
%
% Args:
%
%    db (ts object): time series object
%
%    xrange: date range. Check XXXXX for correct format
%
%    M (integer): number of bins (default 10)
%
%    X (vector): histogram with bin centers given by X
%
% Returns:
%    :
%    h (figure handle): handle to the plotted histogram
%
% Example:
%    ::
%
%       hist('1994m7:1997m1',db(:,:,1),...
%           'figsize',[2,2],...
%           'figtitle','no title',...
%           'logy',true,...
%           'subplots',true)
%
% Note:
%     In addition to those matlab properties, RISE adds further properties,
%     which allow to control for. See parse_plot_args
%

tmp=utils.plot.myplot(@hist,varargin{:});

if nargout
    h=tmp;
end

end
