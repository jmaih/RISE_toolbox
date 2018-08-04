function h=hist(varargin)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


%     h = hist(db) bins the elements of db into 10 equally spaced containers
%     and returns the number of elements in each container.  If db is a
%     matrix, hist works down the columns.
%  
%     h = hist(xrange,db),selects the data in dates xrange
%  
%     h = hist(...,db,M), where M is a scalar, uses M bins.
%  
%     h = hist(...,db,X), where X is a vector, returns the distribution of db
%     among bins with centers specified by X. The first bin includes
%     data between -inf and the first center and the last bin
%     includes data between the last bin and inf. Note: Use HISTC if
%     it is more natural to specify bin edges instead. 
%  
%     In addition to those matlab properties, RISE adds further properties,
%     which allow to control for. See parse_plot_args

% hist('1994m7:1997m1',db(:,:,1),...
%     'figsize',[2,2],...
%     'figtitle','no title',...
%     'logy',true,...
%     'subplots',true)

tmp=utils.plot.myplot(@hist,varargin{:});

if nargout
    h=tmp;
end

end
