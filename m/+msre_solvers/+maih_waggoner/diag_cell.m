function Cout=diag_cell(Cin)
% INTERNAL FUNCTION
%

[h,h2]=size(Cin);

if h~=h2
    
    error('cell array should be square')
    
end

Cout=cell(1,h);

for ireg=1:h
    
    Cout{ireg}=Cin{ireg,ireg};
    
end

end