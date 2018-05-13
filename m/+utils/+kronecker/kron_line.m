function A=kron_line(Acell,order)

% kroneckers of numbers : test=kron_line(1:4,3);
% kroneckers of strings : test=kron_line({'a','b','c','d'},3);

[nr,nc]=size(Acell);

if nr>1 && nc>1
    
    error('this function only takes kroneckers of vectors')
    
end

Acell=Acell(:).';

A=Acell;

cellstr=iscellstr(A);

n=numel(Acell);

for ii=2:order
    
    A=repmat(A(:),1,n);
    
    multiply()
    
    A=A.';
    
    A=A(:).';
    
end

    function multiply()
        
        if cellstr
            
            for jj=1:size(A,1)
                
                A(jj,:)=strcat(A(jj,:),Acell);
                
            end
            
        else
            
            A=bsxfun(@times,A,Acell);
            
        end
        
    end

end