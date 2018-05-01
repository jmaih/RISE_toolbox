function obj=subsasgn(obj,s,b)% subsasgn
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

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
            
            [date_numbers,datta,rows_dates,varloc,pages]=...
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
            
            obj=ts(date_numbers(:),datta,obj.varnames,obj.description);
            
        otherwise
            
            error(['unexpected type "',s.type,'"'])
            
    end
    
end

end
