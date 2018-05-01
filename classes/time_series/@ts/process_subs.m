function [date_numbers,datta,rows_dates,varloc,pages]=process_subs(obj,subs,caller)
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


nsubs=numel(subs);
date_numbers=obj.date_numbers;
datta=obj.data;
varloc=locate_variables(obj.varnames,obj.varnames,true);
rows_dates=1:obj.NumberOfObservations;
pages=1:obj.NumberOfPages;

switch nsubs
    case 1 % variables, dates or time shift
        subs1=subs{1};
        if isnumeric(subs1) % dates or time shift
            freq=serial2frequency(subs1);
            if ~isempty(freq)
                rows_dates=select_from_serial(subs1,freq);
            else
                if strcmp(caller,'subsasgn')
                    error('assignments with shifted dates not allowed')
                end
                lead_or_lag=abs(subs1);
                if subs1<0
                    datta=datta(1:end-lead_or_lag,:,:);
                    date_numbers=date_numbers(1+lead_or_lag);
                elseif subs1>=0
                    datta=datta(1+lead_or_lag:end,:,:);
                    date_numbers=date_numbers(1);
                end
                nobs=size(datta,1);
                rows_dates=1:nobs;
                date_numbers=date_numbers+(0:nobs-1);
            end
        elseif islogical(subs1) % select through the rows
            if (size(subs1,1)~=obj.NumberOfObservations) %&& ~strcmp(caller,'subsasgn')
                error('# elements in the logical call must match # observations')
            end
            rows_dates=subs1;
        elseif ischar(subs1)|| iscellstr(subs1) % dates or variable
            if ischar(subs1)
                subs1=cellstr(subs1); 
            end
            subs1=cellfun(@(x)x(~isspace(x)),subs1,'uniformOutput',false);
            % if starts with a digit: date
            date_style=false;
            if ~isempty(subs1)
                % check that the first element is a letter
                date_style=~isstrprop(subs1{1}(1),'alpha'); % date_style=~isvarname(subs1{1});
            end
            if date_style
                ss=char2serial(subs1);
                freq=serial2frequency(ss); 
                % the call above returns a scalar. so we have to
                % reconstruct a vector for freq in case ss is not a scalar
                freq=freq(1)*ones(size(ss));
                rows_dates=select_from_serial(ss,freq);
            else
                varloc=process_variables(subs1);
            end
        else
            error('wrong format for dates and/or variables')
        end
    case {2,3} % dates and variables
        dates=subs{1};
        if ischar(dates)
            dates(isspace(dates))=[];
        end
        if ischar(dates) && strcmp(dates,':')
            % do nothing as we are taking all the dates
        else
            freq=serial2frequency(dates);
            if isempty(freq)
                dates=char2serial(dates);
                freq=serial2frequency(dates);
            end
            rows_dates=select_from_serial(dates,freq);
        end
        variables=subs{2};
        if isnumeric(variables)
            varloc=variables;
            nvars=obj.NumberOfVariables;
            if ~all(ismember(varloc,1:nvars))
                error(['variables position must be in range [1,',sprintf('%0.0f',nvars),']'])
            end
        else
            varloc=process_variables(variables);
        end
        if nsubs==3 %and pages
            pages_=subs{3};
            if ischar(pages_)
                pages_(isspace(pages_))=[];
                if ~strcmp(pages_,':')
                    error('if pages is a char, then it must be ":"')
                end
            else
                pages=pages_;
            end
        end
    otherwise
        error([caller,' only supports 1 to 3 arguments'])
end

    function varloc=process_variables(vnames)
        done=false;
        if ischar(vnames)
            vnames(isspace(vnames))=[];
            all_vars=strcmp(vnames,':');
            if all_vars
                varloc=1:obj.NumberOfVariables;
                done=true;
            else
                vnames=regexp(vnames,',','split');
            end
        end
        if ~done
            varloc=locate_variables(vnames,obj.varnames);
        end
    end

    function [rows_dates]=select_from_serial(x,freq)
        freq0=serial2frequency(obj.date_numbers(1));
        rows_dates=false(1,obj.NumberOfObservations);
        % do the ones with the same frequency
        %------------------------------------
        good=freq==freq0;
        same_freq=find(good);
        for idate=same_freq
            tmp=find(abs(x(idate)-obj.date_numbers)<1e-9);%<--tmp=find(x(idate)==obj.date_numbers);
            if isempty(tmp)
                error(['date ',char(serial2date(x(idate))),' is out of range'])
            end
            rows_dates(tmp)=true;
        end
        % do the ones with annual frequency
        %----------------------------------
        wrong_freq=find(~good);
        
        if ~isempty(wrong_freq)
            x=x(wrong_freq);
            freq=freq(wrong_freq);
            if ~all(freq==1)
                error('off frequencies can only be annual')
            end
            % it is annual and so get all the years
            year0=floor(obj.date_numbers/freq0);
            year1=floor(x./freq);
            for ii=1:numel(x)
                pointer=year0==year1(ii);
                if ~any(pointer)
                    error(['year "',sprintf('%0.0f',year1(ii)),'" is out of range'])
                end
                rows_dates(pointer)=true;
            end
        end
    end
end
