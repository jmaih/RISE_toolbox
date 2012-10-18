classdef rise_date
    % TODO:
    % 1- define the function end in rise_date
    properties(SetAccess = private, Hidden = true)
        frequency_code
        week_end_days
        last_business_day
    end
    properties
        date
        date_number
        freq
    end
    methods
        function this=rise_date(date,weekend_days,last_business_day)
            % Anuual or yearly: 1900 or '1900'
            % Semi-anuual : '1900H1'
            % Quarterly: '1900Q1'
            % Monthly: '1900M1'
            % Weekly: '1900M1D3W' or '19000103W'
            % Daily: '1900M1D3' or '19000103D'
            % default weekend_days=2. Could also be introduced as a logical vector
            % where each element assumes a value of 1 if the corresponding day
            % is a business day and zero otherwise
            if nargin
                if isa(date,'rise_date')
                    this=date;
                else
                    if nargin<3
                        last_business_day=6;
                    else
                        if ischar(last_business_day)
                            BD={'SUN','MON','TUE','WED','THU','FRI','SAT'
                                1,2,3,4,5,6,7};
                            loc=find(strcmpi(last_business_day,BD(1,:)));
                            if isempty(loc)
                                error([mfilename,':: last business days should be one of ''SUN'',''MON'',''TUE'',''WED'',''THU'',''FRI'',''SAT'''])
                            end
                            last_business_day=BD{2,loc};
                        elseif isnumeric(last_business_day)
                            if ~isscalar(last_business_day)||~ismember(last_business_day,(1:7))
                                error([mfilename,':: last business day expected to be a scalar in [1,7]'])
                            end
                        else
                            error([mfilename,':: last business day must be a scalar in [1,7] or one of ''SUN'',''MON'',''TUE'',''WED'',''THU'',''FRI'',''SAT'''])
                        end
                    end
                    if nargin<2
                        weekend_days=[];
                    end
                    if ~isempty(weekend_days)
                        if isnumeric(weekend_days) && numel(weekend_days)==7 && all(ismember(weekend_days,[0,1]))
                            weekend_days=logical(weekend_days);
                        end
                        if islogical(weekend_days)
                            weekend_days=weekend_days(:)';
                            if numel(weekend_days)~=7
                                error([mfilename,':: weekend_days should be a 7-element vector'])
                            end
                            test=sum(abs(bsxfun(@minus, weekend_days,...
                                [1,0,0,0,0,0,1 % sun & sat
                                1,0,0,0,0,0,0 % sun
                                0,0,0,0,0,0,0 % no wkn
                                0,0,0,0,0,0,1 % sat
                                0,0,0,0,0,1,0])),2); % fri
                            if ~any(test==0) 
                                error([mfilename,':: weekend_days vector not recognized'])
                            end
                        elseif isnumeric(weekend_days)
                            weekend_days=unique(weekend_days);
                            if ~all(ismember(weekend_days,(1:7)))
                                error([mfilename,':: weekend days should be in [1,7]'])
                            end
                            wkn_days=false(1,7);
                            wkn_days(weekend_days)=true;
                            weekend_days=wkn_days;
                        else
                            error([mfilename,':: unrecognized entry for weekend_days'])
                        end
                    end
                    if ~iscellstr(date)
                        if isnumeric(date)
                            date=cellstr(num2str(date(:)));
                        elseif ischar(date)
                            date=cellstr(date);
                        elseif iscell(date)
                            date=cell2mat(date);
                            date=cellstr(num2str(date(:)));
                        end
                    end
                    this=rise_date.empty(0,1);
                    date=upper(date);
                    [this(1).date,this(1).date_number,this(1).freq,this(1).frequency_code]=string_decomposition(date{1});%
                    if ~isempty(weekend_days)&& this(1).frequency_code>4
                        this(1).week_end_days=find(weekend_days);
                        this(1).last_business_day=last_business_day;
                        if ismember(rise_date.weekday(this(1).date_number),this(1).week_end_days)
                            error([mfilename,':: weekly or daily dates occurring on non-business days'])
                        end
                    else
                        if this(1).frequency_code<=4 % annual, semi-annual, quarterly, monthly
                            this(1).week_end_days=find(false(1,7));
                            this(1).last_business_day=1;
                        else % weekly, daily
                            this(1).week_end_days=find(logical([1,0,0,0,0,0,1]));
                            this(1).last_business_day=6;
                        end
                    end
                    is_regular=true;
                    for ii=2:numel(date)
                        this(ii,1)=rise_date(date(ii),this(1).week_end_days,this(1).last_business_day);
                        if this(ii)<=this(ii-1)
                            error([mfilename,':: dataset not monotonically increasing. Break occurs at observation (',int2str(ii),') in ',this(ii).date])
                        end
                        if strcmp(this(1).freq,'W') && ...
                                ~isequal(rise_date.weekday(this(ii).date_number),...
                                rise_date.weekday(this(1).date_number))
                            is_regular=false;
                        end
                    end
                    if ~is_regular
                        warning([mfilename,':: weekly dates are irregular. you may want to consider changing frequency to daily']) %#ok<WNTAG>
                    end
                end
            end
        end
        function [dd,dn]=observation_2_date(this,obs)
            obs=obs(:);
            dn=this.date_number*ones(numel(obs),1);
            positives=obs>0;
            negatives=obs<0;
            % adjustment for negative observations
            obs(negatives)=obs(negatives)-1;
            % now run positives and negatives in turns, not simultaneously
            dn(positives)=sequential_date2obs(dn(positives),obs(positives),+1);
            dn(negatives)=sequential_date2obs(dn(negatives),obs(negatives),-1);
            dv=datevec(dn);
            y=int2str(dv(:,1));
            m=int2str(dv(:,2));
            switch this.freq
                case ''
                    dd=y;
                case 'H'
                    h=int2str(dv(:,2)/6);
                    dd=strcat(y,'H',h);
                case 'Q'
                    q=int2str(dv(:,2)/3);
                    dd=strcat(y,'Q',q);
                case 'M'
                    dd=strcat(y,'M',m);
                case {'W','D'}
                    d=int2str(dv(:,3));
                    dd=strcat(y,'M',m,'D',d);
                    if strcmp(this.freq,'W')
                        dd=strcat(dd,'W');
                    end
            end
            function dnobs=sequential_date2obs(dnobs,obs,sgn)
                iter=1;
                bad=iter<abs(obs);
                while any(bad)
                    iter=iter+1;
                    dnobs(bad)=rise_date.update_date_number(dnobs(bad),this.freq,this.week_end_days,sgn);
                    bad=iter<abs(obs);
                end
            end
        end
        function obs=date_2_observation(this,dd)
            if ~iscell(dd)
                if ischar(dd)
                    dd=cellstr(dd);
                else
                    dd=num2cell(dd);
                end
            end
            obs=nan(numel(dd),1);
            for ii=1:numel(dd)
                sequence=colon(this,dd{ii});
                obs(ii)=numel(sequence);
            end
        end
        function this=colon(this1,this2)
            this2=rise_date(this2);%,weekend_days
            if numel(this2)>1
                error([mfilename,':: dates must be scalar'])
            end
            if this2<this1
                error([mfilename,':: second date cannot occur before first date'])
            end
            if ~isequal(this1.freq,this2.freq)
                error([mfilename,':: inputs must have the same frequency'])
            end
            dn0=this1.date_number;
            dn1=this2.date_number;
            if isequal(this1.frequency_code,5)% strcmp(this1.freq,'W')
                t1=rise_date.weekday(dn0);
                if ~isequal(t1,rise_date.weekday(dn1))
                    error([mfilename,':: weekly dates not observed on the same day'])
                end
            end
            dn=dn0;
            while ~isequal(dn(end),dn1)
                dn=[dn;rise_date.update_date_number(dn(end),this1.freq,this1.week_end_days,+1)];%#ok<AGROW> %,weekend_days
                if dn(end)>dn1
                    error([mfilename,':: incrementing just got pass the target. Something wrong in the incrementation algorithm'])
                end
            end
            dv=datevec(dn);
            y=cellstr(num2str(dv(:,1)));
            m=cellstr(num2str(dv(:,2)));
            d=cellstr(num2str(dv(:,3)));
            switch this1.freq
                case ''
                    this=rise_date(y);
                case 'H'
                    h=cellstr(num2str(dv(:,2)/6));
                    this=rise_date(strcat(y,'H',h));
                case 'Q'
                    q=cellstr(num2str(dv(:,2)/3));
                    this=rise_date(strcat(y,'Q',q));
                case 'M'
                    this=rise_date(strcat(y,'M',m));
                case 'W'
                    this=rise_date(strcat(y,'M',m,'D',d,'W'));
                case 'D'
                    this=rise_date(strcat(y,'M',m,'D',d));
                otherwise
                    error([mfilename,':: unrecognized frequency'])
            end
            
        end
        function this=convert_date(this0,newfreq,new_wknd_days,weekly_day)
            error(nargchk(2,4,nargin,'string'))
            if nargin<4
                weekly_day=6; % friday
                if nargin<3
                    new_wknd_days=false(1,7);
                    new_wknd_days([1,7])=true;
                end
            end
            oldfreq=this0(1).frequency_code;
                freqtable={'','H','Q','M','W','D'
                    1,2,3,4,5,6};
            if ischar(newfreq)
                newfreq_id=find(strcmp(newfreq,freqtable(1,:)));
                if isempty(newfreq_id)
                    error([mfilename,':: unknown frequency ''',newfreq,''])
                end
                newfreq=freqtable{2,newfreq_id};
            end
            if isequal(oldfreq,newfreq)
                this=this0;
            else
                if oldfreq>newfreq
                    error([mfilename,':: date conversion from high to low frequency is undefined'])
                end
                % the date vectors remain unchanged as every date is
                % expressed on an end-of-period basis. e.g. all yearly data
                % have month 12 and day 31. Only the date numbers will have
                % to be adjusted for weekly and possibly daily frequencies.
                dn=vertcat(this0.date_number);
                switch newfreq
                    case 5 %'W': all dates should correspond to the day of the week in which the data are sampled
                        if ~isscalar(weekly_day) && ~ismember(weekly_day,(1:7))
                            error([mfilename,':: last business day must be a scalar integer in [1,7]'])
                        end
                        baddates=rise_date.weekday(dn)~=weekly_day;
                        while any(baddates)
                            dn(baddates)=dn(baddates)-1;
                            baddates=rise_date.weekday(dn)~=weekly_day;
                        end
                    case 6 %'D': all dates should correspond to a business day
                        if islogical(new_wknd_days)
                            if numel(new_wknd_days)~=7
                                error([mfilename,':: weekend_days should be a 7-element vector'])
                            end
                            new_wknd_days=find(new_wknd_days);
                        else
                            new_wknd_days=unique(new_wknd_days);
                            if ~ismember(new_wknd_days,(1:7))
                                error([mfilename,':: weekend days should be between in [1,7]'])
                            end
                        end
                        baddates=ismember(rise_date.weekday(dn),new_wknd_days);
                        while any(baddates)
                            dn(baddates)=dn(baddates)-1;
                            baddates=ismember(rise_date.weekday(dn),new_wknd_days);
                        end
                end
                % penultimate step: date numbers to date vectors
                dato=rise_date.date2date(dn,newfreq);
                this=rise_date(dato);
            end
        end
        function flag=eq(this1,this2)
            this2=rise_date(this2);
            flag=this1.date_number==this2.date_number;
        end
        function flag=ge(this1,this2)
            this2=rise_date(this2);
            flag=this1.date_number>=this2.date_number;
        end
        function flag=gt(this1,this2)
            this2=rise_date(this2);
            flag=this1.date_number>this2.date_number;
        end
        function flag=le(this1,this2)
            this2=rise_date(this2);
            flag=this1.date_number<=this2.date_number;
        end
        function flag=lt(this1,this2)
            this2=rise_date(this2);
            flag=any([this1.date_number]<this2.date_number);
        end
        function this=max(this1,this2)
            this2=rise_date(this2);
            this=this1;
            if this1<this2
                this=this2;
            end
        end
        function this=min(this1,this2)
            this2=rise_date(this2);
            this=this1;
            if this1>this2
                this=this2;
            end
        end
        function this=plus(this1,this2)
            if isa(this2,'double')% && isscalar(this2)
                this=this1.observation_2_date(this2+1);
                this=rise_date(this);
            else
                error([mfilename,':: only defined for adding a rise_date with a double'])
            end
        end
        function this=minus(this1,this2)
            this=plus(this1,-(this2+1));
        end
    end
    methods(Static)
        function dd=weekday(date_number)
            dd=mod(fix(date_number)-2,7)+1;
        end
        function dn=update_date_number(dn0,freq,weekend_days,sgn)
            dn=dn0;
            dv0=datevec(dn0);
            if ~isscalar(sgn)
                error([mfilename,':: sgn argument should be a scalar'])
            end
            increase=sgn>0;
            switch freq
                case {1,''}
                    yy=365;
                    if increase
                        dn=dn0+yy;
                        dv=datevec(dn);
                        leap=rise_date.isleap(dv(:,1));
                        dn(leap)=dn(leap)+1;
                    else
                        dn=dn0-yy;
                        leap=rise_date.isleap(dv0(:,1));
                        dn(leap)=dn(leap)-1;
                    end
                case {2,'H'}
                    hh=[181,184];
                    h1=dv0(:,2)/6==1;
                    h2=dv0(:,2)/6==2;
                    if increase
                        dn(h1,1)=dn0(h1,1)+hh(2);
                        dn(h2,1)=dn0(h2,1)+hh(1);
                        dv=datevec(dn);
                        more=h2 & rise_date.isleap(dv(:,1));
                        dn(more,1)=dn(more,1)+1;
                    else
                        dn(h1,1)=dn0(h1,1)-hh(1);
                        dn(h2,1)=dn0(h2,1)-hh(2);
                        more=h1 & rise_date.isleap(dv0(:,1));
                        dn(more,1)=dn(more,1)-1;
                    end
                case {3,'Q'}
                    qq=[90,91,92,92];
                    dv=datevec(dn0);
                    q1=dv(:,2)/3==1;
                    q2=dv(:,2)/3==2;
                    q3=dv(:,2)/3==3;
                    q4=dv(:,2)/3==4;
                    if increase
                        dn(q1,1)=dn0(q1,1)+qq(2);
                        dn(q2|q3,1)=dn0(q2|q3,1)+qq(4);
                        dn(q4,1)=dn0(q4,1)+qq(1);
                        dv=datevec(dn);
                        more=q4 & rise_date.isleap(dv(:,1));
                        dn(more,1)=dn(more,1)+1;
                    else
                        dn(q2,1)=dn0(q2,1)-qq(2);
                        dn(q3|q4,1)=dn0(q3|q4,1)-qq(4);
                        dn(q1,1)=dn0(q1,1)-qq(1);
                        more=q1 & rise_date.isleap(dv0(:,1));
                        dn(more,1)=dn(more,1)-1;
                    end
                case {4,'M'}
                    mm=[31,28,31,30,31,30,31,31,30,31,30,31];
                    m1=dv0(:,2)==1; m2=dv0(:,2)==2; m3=dv0(:,2)==3;
                    m4=dv0(:,2)==4; m5=dv0(:,2)==5; m6=dv0(:,2)==6;
                    m7=dv0(:,2)==7; m8=dv0(:,2)==8; m9=dv0(:,2)==9;
                    m10=dv0(:,2)==10; m11=dv0(:,2)==11; m12=dv0(:,2)==12;
                    if increase
                        m30=m3|m5|m8|m10;
                        dn(m30,1)=dn0(m30,1)+mm(4);
                        m31=m2|m4|m6|m7|m9|m11|m12;
                        dn(m31,1)=dn0(m31,1)+mm(1);
                        dn(m1,1)=dn0(m1,1)+mm(2);
                        more=m1 & rise_date.isleap(dv0(:,1));
                        dn(more,1)=dn(more,1)+1;
                    else
                        m30=m4|m6|m9|m11;
                        dn(m30,1)=dn0(m30,1)-mm(4);
                        m31=m1|m3|m5|m7|m8|m10|m12;
                        dn(m31,1)=dn0(m31,1)-mm(1);
                        dn(m2,1)=dn0(m2,1)-mm(2);
                        more=m2 & rise_date.isleap(dv0(:,1));
                        dn(more,1)=dn(more,1)-1;
                    end
                case {5,'W'}
                    ww=7;
                    if increase
                        dn=dn0+ww;
                    else
                        dn=dn0-ww;
                    end
                case {6,'D'}
                    dd=1;
                    if increase
                        dn=dn0+dd;
                        bad=ismember(rise_date.weekday(dn),weekend_days);
                        while any(bad)
                            dn(bad)=dn(bad)+dd;
                            bad=ismember(rise_date.weekday(dn),weekend_days);
                        end
                    else
                        dn=dn0-dd;
                        bad=ismember(rise_date.weekday(dn),weekend_days);
                        while any(bad)
                            dn(bad)=dn(bad)-dd;
                            bad=ismember(rise_date.weekday(dn),weekend_days);
                        end
                    end
                otherwise
                    error([mfilename,':: unrecognized frequency ''',freq,''])
            end
        end
        function flag=isleap(x)
            flag=(mod(x,4)==0 & mod(x,100)~=0)|...
                (mod(x,4)==0 & mod(x,100)==0 & mod(x,400)==0);
        end
        function date=date2date(datein,frequency_code)
            is_date_vector=size(datein,2)==3;
            if ~is_date_vector
                datein=datevec(datein);
            end
            y=datein(:,1);
            m=datein(:,2);
            d=datein(:,3);
            switch frequency_code
                case 1% if is_yearly
                    date=int2str(y);
                case 2% elseif is_semiannual
                    date=strcat(int2str(y),'H',int2str(m/6));
                case 3% elseif is_quarterly
                    date=strcat(int2str(y),'Q',int2str(m/3));
                case 4% elseif is_monthly
                    date=strcat(int2str(y),'M',int2str(m));
                case {5,6}% elseif is_weekly || is_daily
                    date=strcat(int2str(y),'M',int2str(m),'D',int2str(d));
                    if frequency_code==5% is_weekly
                        date=strcat(date,'W');
                    end
            end
        end
    end
end
function [date,date_nbr,freq,frequency_code]=string_decomposition(string)
string(isspace(string))=[];
M=strfind(string,'M');
H=strfind(string,'H');
Q=strfind(string,'Q');
W=strfind(string,'W');
D=strfind(string,'D');
monthsdays=[31,28,31,30,31,30,31,31,30,31,30,31];
%           J   F  M  A  M  J  J  A  S  O  N  D

y=str2double(string);
if ~isnan(y)
    if ~isequal(y,floor(y))
        error([mfilename,':: wrong date format ''',string,''])
    end
    m=12;
    d=monthsdays(m); % december is a 31-day month
    frequency_code=1;% is_yearly=true;
    freq='';
else
    if ~isempty(H)
        y=str2double(string(1:H-1));
        m=6*str2double(string(H+1:end));
        if isnan(y)|| ~isequal(y,floor(y))|| ~ismember(m,(1:12))
            error([mfilename,':: wrong semi-annual date format ''',string,''])
        end
        d=monthsdays(m);
        frequency_code=2;% is_semiannual=true;
        freq='H';
    else
        if ~isempty(Q)
            y=str2double(string(1:Q-1));
            m=3*str2double(string(Q+1:end));
            if isnan(y)|| ~isequal(y,floor(y))||~ismember(m,(1:12))
                error([mfilename,':: wrong quarterly date format ''',string,''])
            end
            d=monthsdays(m);
            frequency_code=3;% is_quarterly=true;
            freq='Q';
        else
            if ~isempty(M)
                y=str2double(string(1:M-1));
                if isnan(y)|| ~isequal(y,floor(y))
                    error([mfilename,':: wrong monthly, weekly or daily date format ''',string,''])
                end
                if rise_date.isleap(y)
                    monthsdays(2)=29;
                end
                m=str2double(string(M+1:end));
                if ~isnan(m)
                    d=monthsdays(m);
                    frequency_code=4;% is_monthly=true;
                    freq='M';
                else
                    if ~isempty(D)
                        m=str2double(string(M+1:D-1));
                        if isnan(m)||~ismember(m,(1:12))
                            error([mfilename,':: wrong monthly, weekly or daily date format ''',string,''])
                        end
                        d=str2double(string(D+1:end));
                        if isnan(d)
                            if strcmp(string(end),'W')
                                d=str2double(string(D+1:W-1));
                                if isnan(d)|| ~ismember(d,(1:monthsdays(m)))
                                    error([mfilename,':: wrong weekly date format ''',string,''])
                                end
                                frequency_code=5;% is_weekly=true;
                                freq='W';
                            else
                                error([mfilename,':: wrong weekly or daily date format ''',string,''])
                            end
                        else
                            if ~ismember(d,(1:monthsdays(m)))
                                error([mfilename,':: wrong daily date format ''',string,''])
                            end
                            frequency_code=6;% is_daily=true;
                            freq='D';
                        end
                    else
                        error([mfilename,':: wrong monthly, weekly or daily date format ''',string,''])
                    end
                end
            else
                if length(string)<6
                    error([mfilename,':: wrong weekly or daily date format ''',string,''])
                end
                d=str2double(string(end-2:end-1));
                m=str2double(string(end-4:end-3));
                y=str2double(string(1:end-5));
                if rise_date.isleap(y)
                    monthsdays(2)=29;
                end
                if any(isnan([y,m,d]))|| ~isequal(y,floor(y))|| ~ismember(m,(1:12))|| ~ismember(d,(1:monthsdays(m)))
                    disp(help('rise_date'))
                    error([mfilename,':: wrong weekly or daily date format ''',string,''''])
                end
                if strcmp(string(end),'D')
                    frequency_code=6;% is_daily=true;
                    freq='D';
                elseif strcmp(string(end),'W')
                    frequency_code=5;% is_weekly=true;
                    freq='W';
                else
                    error([mfilename,':: wrong weekly or daily date format ''',string,''])
                end
            end
        end
    end
end
date_nbr=datenum([y,m,d]);
if size(date_nbr,2)>1
    disp([y,m,d])
    error([mfilename,':: date number calculation failed. Your dates may be too large...'])
end

date=rise_date.date2date([y,m,d],frequency_code);

end