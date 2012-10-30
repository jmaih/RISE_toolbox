function flag=is_atom(string)
flag=isempty(regexp(string,'[/*\-+^()]','start'));
end

