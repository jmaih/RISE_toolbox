function [s,frequency]=date2serial(dat,silent)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

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

frequency=[];

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

if ischar(dat)
    % in order to take the appropriate size below
    dat=cellstr(dat);
    
end

s=reshape(s,size(dat));


    function do_it_then()
        
        [dec,flag]=decompose_date(dat);
        
        if ~flag
            
            error('decomposition failed')
            
        end
        
        year=cellfun(@str2double,{dec.year},'uniformOutput',true);
        
        period=cellfun(@str2double,{dec.period},'uniformOutput',true);
        
        frequency={dec.WMQH};
        
        freq=frequency2num(frequency);
        
        stamp=time_frequency_stamp();
        
        dn=@(x,per,freq)x.*freq+per-1+stamp(freq);
        
        s=dn(year,period,freq);
        
    end


end

