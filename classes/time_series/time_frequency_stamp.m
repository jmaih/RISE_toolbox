function [stamp,unstamp]=time_frequency_stamp(option)
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

% stamps helps put a stamp on the serial numbers such that the frequency is
% recoverable

% unstamp inverts the stamp above

if nargin==0
    
    option=2;
    
elseif ~ismember(option,[1,2])
    
    error('case not implemented')
    
end

stamp=@mystamp;

unstamp=@myunstamp;

    function us=myunstamp(x)
        
        if option==1
            
            us=round(100*x);
            
        elseif option==2
            
            us=round(1./x-1);
                        
        end
        
    end

    function s=mystamp(x)
        
        if option==1
            
            s=0.01*x;
            
        elseif option==2
            
            s=1./(1+x);
            
        end
        
    end

end