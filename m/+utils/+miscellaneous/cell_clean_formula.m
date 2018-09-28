function C=cell_clean_formula(C0)

% (X)^2--> X^2
C1=regexprep(C0,'\((\w+)\)\^2','$1^2');

%(X)*(Y)--> X*Y
C2=regexprep(C1,'\((\w+)\)\*\((\w+)\)','$1*$2');

% -(X)--> -X
C3=regexprep(C2,'\((-)?([\w+\.]+)\)','$1$2');

% c*(A+B) --> c*A+c*B
myreplace=@another_slick; %#ok<NASGU>

C4=regexprep(C3,'(-)?(\d+[\d+.]*)\*\(([^\(\)]+)\)','${myreplace($1,$2,$3)}');
% C4=regexprep(C3,'(-)?([\d+\.]+)\*\(([^\(\)]+)\)','${myreplace($1,$2,$3)}');

C5=strrep(C4,'+-','-');

% (T^2)-->T^2
C=regexprep(C5,'\((\w+\^\d+)\)','$1');

    function d=another_slick(a,b,c)
                
        ab=[a,b];
        
        if strcmp(b,'0')
            
            d=ab;
            
            return
            
        end
        
        abneg=ab(1)=='-';
        
        if abneg
            
            ab=ab(2:end);
            
        end
        
        first_neg=c(1)=='-';
        
        if ~first_neg
            
            c=['+',c];
            
        end
        
        d=regexprep(c,'(\+|\-)',['$1',ab,'*']);
        
        if abneg % flip signs
            
            d=strrep(d,'+','#');
            
            d=strrep(d,'-','@');
            
            d=strrep(d,'#','-');
            
            d=strrep(d,'@','+');
            
        end
        
        if d(1)=='+'
            
            d=d(2:end);
            
        end
        
    end

end