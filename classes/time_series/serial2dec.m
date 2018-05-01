function [dec]=serial2dec(x,checked)
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

if nargin<2
    
    checked=false;
    
end

if ~checked && ~is_serial(x)
    
    error('input must be serial')
    
end

nx=numel(x);

[stamp,unstamp]=time_frequency_stamp();

freq=unstamp(x-floor(x));

year=floor(x./freq);

period=round(x-year.*freq+1-stamp(freq));

dec=initialize();

for ii=1:nx
    
    dec(ii).year=year(ii);
    
    dec(ii).period=period(ii);
    
end

% dat=year_str;
% 
% period_str=cellstr(int2str(period(:)));
% 
% year_str=cellstr(int2str(year(:)));
% 
% frequency=cellstr(frequency2char(freq));


    function dec=initialize()
        
        dec=decompose_date();
        
        dec.freq=freq(1);
        
        dec.frequency=frequency2char(freq(1));
        
        dec.period=1;
        
        dec=dec(1,ones(nx,1));
    end

end