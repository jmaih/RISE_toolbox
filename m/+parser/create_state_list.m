function x=create_state_list(name,number)
% INTERNAL FUNCTION
%

if isscalar(number)
    
    number=(1:number).';
    
else
    
    number=number(:);
    
end

x=cellstr(strcat(name,'_',num2str(number)));

x=cellfun(@(z)z(~isspace(z)),x,'uniformOutput',false);

x=x(:).';

x=x(1:number);

end
