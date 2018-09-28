function C=cellmult(A,B,precision)
% Symbolic Multiplication of (matrix OR cellstr) with (matrix OR cellstr)
% designed for use with Ben Petschel stuff on Grobner basis

if nargin<3
    
    precision=20;
    
end

[A,ra,ca]=utils.miscellaneous.cell_format_input(A,precision);

[B,rb,cb]=utils.miscellaneous.cell_format_input(B,precision);

if ca~=rb
    
    if ca==1 && ra==1
        
        C=strcat('(',A{1},')*','(',B,')');
        
    else
        
        error('wrong size for multiplication')
        
    end
    
else
    
    C=cell(ra,cb);
    
    for irow=1:ra
        
        for icol=1:cb
            
            arow=A(irow,:); arow=strcat('(',arow(:),')');
            
            bcol=strcat('(',B(:,icol),')');
            
            same=strcmp(arow,bcol);
            
            tmp=strcat(arow,'*',bcol);
            
            tmp(same)=strcat(arow(same),'^2');
            
            tmp=cell2mat(strcat(tmp(:).','+'));
            
            C{irow,icol}=tmp(1:end-1);
            
        end
        
    end
    
end

C=utils.miscellaneous.cell_clean_formula(C);

end