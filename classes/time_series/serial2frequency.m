function [freq,frequency,year,period]=serial2frequency(x)
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

freq=[];

frequency=[];

year=[];

period=[];

if isnumeric(x)
    
    [stamp,unstamp]=time_frequency_stamp();
    
    freq_=unstamp(x-floor(x));
    
    if ~all(freq_==freq_(1))
        
        error('dates should have the same frequency')
        
    end
    
    fmap=frequency_map;
    
    if any(freq_(1)==fmap.code);
        
        freq=freq_;
        
        if nargout>1
            
            frequency=frequency2char(freq(1));
            
            if nargout>2
                
                year=floor(x./freq);
                
                period=round(x-year.*freq+1-stamp(freq));
                
            end
            
        end
        
    end
    
end

end