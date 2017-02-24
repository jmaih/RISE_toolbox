function opt=disp_defaults(defs)

% remove comments that call for further action:
% e.g. 'irf_shock_list(sr)' into 'irf_shock_list'
%-------------------------------------------------
strings=regexprep(defs(:,1),'\<(\w+)\>.*','$1'); % (\<\w+\>).* would also work

values=defs(:,2);

myopt=cell2struct(values,strings,1);

if nargout
    
    opt=myopt;
    
else
    
    disp(' ')
    
    disp(myopt)
    
    disp(' ')
    
end

end