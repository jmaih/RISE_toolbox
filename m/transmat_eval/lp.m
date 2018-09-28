function v=lp(m,c)

if max(c)>numel(m)
    
    v=m(1);
    
elseif numel(c)==1
    
    v=m(c);
    
else
    
    v=mean(m(c));
    
end

end