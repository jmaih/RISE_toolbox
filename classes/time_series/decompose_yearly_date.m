function [dec,flag]=decompose_yearly_date(x)

if ~iscellstr(x)
    
    error('input must be a cellstr')
    
end

dec=[];

flag=false;

nx=numel(x);

ii=0;

tmp=decompose_date();

tmp=tmp(1,ones(nx,1));

while ii<nx
    
    ii=ii+1;
    
    xx=str2double(x{ii});
    
    flag_i=xx==floor(xx);
    
    if ~flag_i
        
        break
        
    end
    
    tmp(ii).year=x{ii}; tmp(ii).period='1';
    
end

if flag_i
    
    flag=true;
    
    dec=reshape(tmp,size(x));
    
end

end