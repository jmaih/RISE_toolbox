function k=freq2freq(freq)

k=1;

if isempty(freq),return,end

switch freq
    
    case 'Q'
        
        k=4;
        
    case 'H'
        
        k=2;
        
    case 'M'
        
        k=12;
        
    case 'W'
        
        k=52;
        
end

end