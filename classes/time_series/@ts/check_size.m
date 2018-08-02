function s=check_size(data)
% INTERNAL FUNCTION: Confirm that the data is consistent with the number of observed dates and variables
%

sz=size(data);

if numel(sz)>3

    error('ts does not handle time series with more than 3 dimensions at the moment')

end

if nargout

    s=sz;

end

end
