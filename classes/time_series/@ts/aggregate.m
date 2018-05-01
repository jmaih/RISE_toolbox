function this=aggregate(this,newfreq,method)
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

if nargin<3
    
    method='distribution';
    
end

FREQ=this.frequency;

newfreq=upper(newfreq);

if ~ismember(newfreq,{'','Q','M','H','W','D'})
    
    error([mfilename,':: unknown type of frequency ',newfreq])
    
end

if strcmp(FREQ,newfreq)
    
    return
    
end

[first,last,number]=first_last_number(FREQ,newfreq);

the_nans=isnan(this);

if size(the_nans,3)>1
    
    error('cannot handle multiple pages')
    
end

good_rows=~any(the_nans,2);

first_guy=find(good_rows,1,'first');

last_guy=find(good_rows,1,'last');

if isempty(first_guy)||isempty(last_guy)||~all(good_rows(first_guy:last_guy))
    
    error('holes in the data, cannot aggregate')
    
end

% now find a location where the first date matches the first date of the
% new frequency
%-----------------------------------------------------------------------
date_numbers=this.date_numbers;

[year,period]=date2year_period(date_numbers);

nobs=numel(date_numbers);

while ~matches_first(first_guy)
    
    first_guy=first_guy+1;
    
    if first_guy>=nobs
        
        error('insufficient number of observations for aggregation')
        
    end
    
end

% now find a location where the last date matches the last date of the
% new frequency
%-----------------------------------------------------------------------
while ~matches_last(last_guy)
    
    last_guy=last_guy-1;
    
    if last_guy<=1
        
        error('insufficient number of observations for aggregation')
        
    end
    
end

datta=this.data;

switch method
    
    case 'interpolation'
        
        indexfun=@(x)(x-1)*number+1;
        
    case 'distribution'
        
        indexfun=@(x)(x-1)*number+1:x*number;
        
    otherwise
        
        error([mfilename,':: unknown aggregation method ',method])
        
end

T=last_guy-first_guy+1;

q=T/number;

C=zeros(q,T);

coef=1/number;

for ii=1:q
    
    C(ii,indexfun(ii))=coef;
    
end

Cdatta=C*datta(first_guy:last_guy,:);

startdate=year(first_guy);

if ~isempty(newfreq)
    
    period_start=find(period(first_guy)==first);
    
    startdate=sprintf('%0.0f%s%0.0f',startdate,newfreq,period_start);
    
end

this=ts(startdate,Cdatta,this.varnames);

    function flag=matches_first(guy)
        
        flag=period(guy)==first;
        
    end

    function flag=matches_last(guy)
        
        flag=period(guy)==last;
        
    end

end

function [first,last,number]=first_last_number(oldfreq,newfreq)

newfreq=frequency2num(newfreq);

switch oldfreq
    
    case {'M',12}
        
        if newfreq>frequency2num('M')
            
            error([mfilename,':: aggregation to higher frequency impossible'])
            
        end
        
        if newfreq==4 %'Q'
            
            first=[1,4,7,10]; last=[3,6,9,12];
            
        elseif newfreq==2%,'H'
            
            first=[1,7]; last=[6,12];
            
        elseif newfreq==1% '' annual
            
            first=1; last=12;
            
        end
        
    case {'Q',4}
        
        if newfreq>frequency2num('Q')
            
            error([mfilename,':: aggregation to higher frequency impossible'])
            
        end
        
        if newfreq==2%,'H'
            
            first=[1,3]; last=[2,4];
            
        elseif newfreq==1% '' annual
            
            first=1; last=4;
            
        end
        
    case {'H',2}
        
        if newfreq>frequency2num('H')
            
            error([mfilename,':: aggregation to higher frequency impossible'])
            
        end
        
        first=1; last=2;
        
    case {'',1}
        
        if newfreq>frequency2num('')
            
            error([mfilename,':: aggregation to higher frequency impossible'])
            
        end
        
end

number=last(1)-first(1)+1;

end
