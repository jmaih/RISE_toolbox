function C=cellplus(A,B,precision)
% Symbolic Multiplication of (matrix OR cellstr) with (matrix OR cellstr)
% designed for use with Ben Petschel stuff on Grobner basis

if nargin<3
    
    precision=20;
    
end

[A,ra,ca]=utils.miscellaneous.cell_format_input(A,precision);

[B,rb,cb]=utils.miscellaneous.cell_format_input(B,precision);

if ca==1 && ra==1
    
    A=repmat(A,size(B));
    
elseif rb==1 && cb==1
    
    B=repmat(B,size(A));
    
end
    
C=strcat(A,'+',B);

C=utils.miscellaneous.cell_clean_formula(C);

end