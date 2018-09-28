function O=cell_str2cell(str)

% Utility function that transforms a string whose evaluation would give a
% matrix, into a cell array.
%
% designed for use with Ben Petschel stuff on Grobner basis

str=strrep(str,'[','');

str=strrep(str,']','');

O=str;

for ii=1:numel(str)
    
    stri=str{ii};
    
    mysplit=regexp(stri,';','split');
    
    nrows=numel(mysplit);
    
    ncols=sum(mysplit{1}==',')+1;
    
    Oi=cell(nrows,ncols);
    
    for irow=1:nrows
        
        Oi(irow,:)=regexp(mysplit{irow},',','split');
        
    end
    
    O{ii}=Oi;
    
end

end