function x=cellstringize(x)

if ischar(x)
    
    x=cellstr(x);
    
end

if ~iscellstr(x)
    
    error('names must be a char or a cellstr')
    
end

good=cellfun(@isvarname,x);

badx=x(~good);

if ~isempty(badx)
    
    disp(badx)
    
    error('the above variables are not valid variable names')
    
end

end
