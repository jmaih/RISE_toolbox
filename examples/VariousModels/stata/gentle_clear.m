function gentle_clear()

noclearall=evalin('base','exist(''noclearall'',''var'')');

if noclearall
    
    noclearall=evalin('base','noclearall');
    
else
    
    noclearall=false;
    
end

if ~noclearall

evalin('base','clear')

end

end