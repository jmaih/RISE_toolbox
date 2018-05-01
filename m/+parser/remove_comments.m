function rawline_=remove_comments(rawline_)
% remove_comments - removes comments of type "//", "%"
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

treat_double_quotes()

% double slash
%--------------
loc_=strfind(rawline_,'//');

if ~isempty(loc_)
    
    rawline_=rawline_(1:loc_(1)-1);
    
end

% percent
%--------
loc_=strfind(rawline_,'%');

if ~isempty(loc_)
    
    rawline_=rawline_(1:loc_(1)-1);
    
end

undo_double_quotes()

    function undo_double_quotes()
        
        rawline_=strrep(rawline_,'¤','%');
        
        rawline_=strrep(rawline_,'§','/');
        
    end

    function treat_double_quotes()
        
        dbl_quotes=strfind(rawline_,'"');
        
        while ~isempty(dbl_quotes)
            
            if numel(dbl_quotes)==1
                
                break
                
            end
            
            quotes=dbl_quotes(1:2);
            
            stretch=quotes(1):quotes(2);
            
            % should replace with exactly one character
            rawline_(stretch)=strrep(rawline_(stretch),'%','¤');
            
            % should replace with exactly two character
            rawline_(stretch)=strrep(rawline_(stretch),'//','§§');
            
            dbl_quotes=dbl_quotes(3:end);
            
        end
        
    end

end
