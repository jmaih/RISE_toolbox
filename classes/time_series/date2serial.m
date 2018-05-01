function [s,frequency]=date2serial(dat,silent)
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

% test=(date2serial(1990):date2serial(1995))
% serial2date(test)
% test=(date2serial('1990'):date2serial('1995'))
% serial2date(test)
% test=(date2serial('1990h1'):date2serial('1995h2'))
% serial2date(test)
% test=(date2serial('1990q1'):date2serial('1995q4'))
% serial2date(test)
% test=(date2serial('1990m1'):date2serial('1995m6'))
% serial2date(test)

if nargin<2
    
    silent=false;
    
end

s=[];

frequency='';

if isempty(dat)
    
    return
    
end

% quick exit if it is already serial
%-----------------------------------

if is_serial(dat)
    
    s=dat;
    
    [~,frequency]=serial2frequency(dat);
    
else
    
    do_it_then()
    
end

if ~silent
    
    if iscellstr(frequency)
        
        if ~all(strcmp(frequency,frequency{1}))
            
            error('frequencies are not the same')
            
        end
        
        frequency=frequency{1};
        
    end
    
end

    function do_it_then()
        
        % try turning into annual if possible and do not scream
        %-------------------------------------------------------
        [dat,is_colon]=char2num(dat,true);
        
        if ischar(dat)||iscellstr(dat)
            
            [dec,flag]=decompose_wmqh_date(dat);
            
            if ~flag
                
                error('decomposition failed')
                
            end
            
            frequency=dec(1).frequency;
            
            year=[dec.year];
            
            period=[dec.period];
            
            freq=[dec.freq];
            
        else
            
            if ~(isnumeric(dat) && all(abs(dat-floor(dat))<1e-12))
                
                error('wrong date format')
                
            end
            
            year=dat;
            
            period=1;
            
            freq=1;
            
        end
        
        s=dec2serial(year,period,freq);
        
        if is_colon
            
            s=s(1):1:s(end);
            
        end
        
        s=s(:).';
        
    end


end

