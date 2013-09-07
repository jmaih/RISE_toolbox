function RawFile=read_file(FileName)

RawFile=cell(0,1);
if isempty(FileName)
    return
end
fid = fopen(FileName);
iter=0;
while 1
    rawline = fgetl(fid);
    if ~ischar(rawline), break, end
    rawline=parser.remove_comments(rawline);
    iter=iter+1;
    if all(isspace(rawline))
        continue
    end
    rawline={rawline,FileName,iter}; %#ok<*AGROW>
    RawFile=[RawFile;rawline];
end
fclose(fid);

end
