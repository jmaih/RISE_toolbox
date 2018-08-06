function obj=subsasgn(obj,s,b) 
% subsasgn Subscripted assignment for time series objects
%
% ::
%
%   db(I)= B assigns the values of B into a time series formed from the
%    elements of db specified by the subscript vector I. 
%
%   db(I,J)=B assigns the values of B into a time series formed from the
%    elements of the rectangular submatrix of db specified by the subscript
%    vectors I and J  
%
%   db(I,J,K)=B assigns the values of B into a time series formed from the
%    elements of the cubic array of db specified by the subscript vectors I,
%    J and K  
%
% Args:
%
%    - one index: variables, dates or time shift
%
%    - two indexes: dates (first dimension) and variables (second
%       dimension)
%
%    - three indexes: dates (first dimension), variables (second dimension)
%       and pages (third dimension) 
%
% Returns:
%    :
%
%    - **this** [ts] : time series
%
% More info:
%    :
%
%    - **specification of dates** : 
%       - list : '1991Q3', '1990Q1:1990Q3,1991,1992Q2', etc. 
%       - serial numbers : date2serial('2011Q1'):date2serial('2014Q1')
%       - boolean indexing [deprecated?]: e.g. db{db{'v1'}>0}
%
%    - **specification of variables** : 
%       - list : 'v1', 'v1,v4', {'v1','v4'} 
%       - numbers : 1, [1,4], etc. but only when dates have been specified 
%
%    - **specification of pages** : 
%       - numbers : 1, 2:4, [1,3,5], etc. 
%       - function handles : e.g. @isreal : only pages where all elements
%          are real
%
%    - **specification of time shift** : Not supported
%
% See also ts.subsref

if numel(s)>1 
    
    if ~strcmp(s(1).type,'.')
       
        error('Unknown form of subsasgn')
        
    end
    
    obj.(s(1).subs)=builtin(mfilename,obj.(s(1).subs),s(2:end),b);
    
else
    
    switch s.type
        
        case '.'
            
            % subs = character string
            obj=builtin(mfilename,obj,s,b);
            
        case {'()','{}'}
            
            [~,datta,rows_dates,varloc,pages]=...
                process_subs(obj,s.subs,mfilename);
            
            if isvector(rows_dates)
                
                if isscalar(b)
                    
                    datta(rows_dates,varloc,pages)=b;
                    
                else
                    
                    datta(rows_dates,varloc,pages)=b;
                    
                end
                
            else
                
                datta(rows_dates)=b;
                
            end
            
            obj=set(obj,'data',datta);
            
        otherwise
            
            error(['unexpected type "',s.type,'"'])
            
    end
    
end

end
