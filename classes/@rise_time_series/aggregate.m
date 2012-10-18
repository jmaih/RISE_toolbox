function this=aggregate(this,newfreq,method)
if nargin<3
    method='distribution';
end
FREQ=this.frequency;
newfreq=upper(newfreq);
if ~ismember(newfreq,{'','Q','M','H','W','D'})
    error([mfilename,':: uknown type of frequency ',newfreq])
end
[first,last,number]=first_last_number();

the_dates={this.TimeInfo.date};
first_not_found=true;
last_not_found=true;
iter=0;
datta=double(this);
numobs=this.NumberOfObservations+1;
while first_not_found
    iter=iter+1;
    loc=strfind(the_dates{iter},FREQ);
    period=str2double(the_dates{iter}(loc+1:end));
    if ismember(period,first) && ~any(isnan(datta(iter,:)))
        period_start=int2str(find(period==first));
        year_start=the_dates{iter}(1:loc-1);
        first_=iter;
        first_not_found=false;
    end
    if iter>=numobs
        error([mfilename,':: insufficient number of valid observations'])
    end
end
iter=this.NumberOfObservations+1;
while last_not_found
    iter=iter-1;
    loc=strfind(the_dates{iter},FREQ);
    period=str2double(the_dates{iter}(loc+1:end));
    if ismember(period,last) && ~any(isnan(datta(iter,:)))
        last_=iter;
        last_not_found=false;
    end
    if iter<number
        error([mfilename,':: insufficient number of valid observations'])
    end
end
switch method
    case 'interpolation'
        indexfun=@(x)(x-1)*number+1;
    case 'distribution'
        indexfun=@(x)(x-1)*number+1:x*number;
    otherwise
        error([mfilename,':: unknown aggregation method ',method])
end
T=last_-first_+1;
q=T/number;
C=zeros(q,T);
coef=1/number;
for ii=1:q
    C(ii,indexfun(ii))=coef;
end
Cdatta=C*datta(first_:last_,:);
startdate=year_start;
if ~isempty(newfreq)
    startdate=[startdate,newfreq,period_start];
end
this=rise_time_series(startdate,Cdatta,this.varnames);
    function [first,last,number]=first_last_number()
        switch FREQ
            case 'M'
                if ismember(newfreq,{'W','D'})
                    error([mfilename,':: aggregation to higher frequency impossible'])
                end
                if strcmp(newfreq,'Q')
                    first=[1,4,7,20]; last=[3,6,9,12];
                elseif strcmp(newfreq,'H')
                    first=[1,7]; last=[6,12];
                elseif strcmp(newfreq,'')
                    first=1; last=12;
                end
            case 'Q'
                if ismember(newfreq,{'W','D','M'})
                    error([mfilename,':: aggregation to higher frequency impossible'])
                end
                if strcmp(newfreq,'H')
                    first=[1,3]; last=[2,4];
                elseif strcmp(newfreq,'')
                    first=1; last=4;
                end
            case 'H'
                if ismember(newfreq,{'W','D','M','Q'})
                    error([mfilename,':: aggregation to higher frequency impossible'])
                end
                first=1; last=2;
            case ''
                if ismember(newfreq,{'W','D','M','Q','H'})
                    error([mfilename,':: aggregation to higher frequency impossible'])
                end
        end
        number=last(1)-first(1)+1;
    end
end
