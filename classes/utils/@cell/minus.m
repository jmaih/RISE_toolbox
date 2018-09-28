function c=minus(a,b)
% minus -- overload minus for cell
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

if ischar(a)
    
    a=cellstr(a);
    
end

if ischar(b)
    
    b=cellstr(b);
    
end

if isempty(b)
    
    c=a;
    
elseif isempty(a)
    % same as above: elements of a that are not in b
    
    c=a;
    
else
    
    [c,ia]=setdiff(a,b);
    
    [~,isort]=sort(ia);
    
    c=c(isort);
    
end

end
