function write2file(C,filename)

fid = fopen(filename,'w+');

if fid == -1
    
    error(['Cannot open file ''%s'' for writing.',filename]);
    
end

if iscellstr(C)
    
    C = sprintf('%s\n',C{:});
    
    if ~isempty(C)
        
        C(end) = '';
        
    end
    
end

count = fwrite(fid,C,'char');

if count ~= length(C)
    
    fclose(fid);
    
    error(['Cannot write character string to file ''%s''.',filename]);
    
end

fclose(fid);

end
