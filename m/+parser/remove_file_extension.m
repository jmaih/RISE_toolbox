function fname=remove_file_extension(fname)
thedot=find(fname=='.');
if ~isempty(thedot)
    fname=fname(1:thedot-1);
end
end