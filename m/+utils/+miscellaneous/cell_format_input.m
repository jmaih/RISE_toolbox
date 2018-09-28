function [A,ra,ca]=cell_format_input(A,precision)

[ra,ca]=size(A);

if isa(A,'double')
    
    A=num2cell(A);
    
    A=cellfun(@(x)num2str(x,precision),A,'uniformOutput',false);
    
end

if ~iscellstr(A) %#ok<ISCLSTR>
    
    error('input of wrong format')
    
end

end