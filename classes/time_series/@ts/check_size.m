function s=check_size(data)

sz=size(data);

if numel(sz)>3
    
    error('ts does not handle time series with more than 3 dimensions at the moment')

end

if nargout
    
    s=sz;
    
end

end
