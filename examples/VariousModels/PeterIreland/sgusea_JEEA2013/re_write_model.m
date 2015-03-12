function re_write_model(m,filname)

% this function rewrites the equation of the model in a file
% so that the complicated equations can be inspected more easily

eqtns=regexprep(m.equations.dynamic,'(\{)(\d+)(\})','$1+$2$3');

fid=fopen([filname,'.rs'],'w');

for ieqtn=1:numel(eqtns)
    eqtn=eqtns{ieqtn};
    fprintf(fid,'\n\n %s',eqtn);
end

fclose(fid);

end