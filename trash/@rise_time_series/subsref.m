function this=subsref(obj,s)
this=obj;
while ~isempty(s)
    if isequal(class(this),'rise_time_series')
        if strcmp(s(1).type,'.')
            this=builtin('subsref',this,s);
            s=[];
            %                     elseif  strcmp(s(1).type,'()')
            %                         if ischar(s(1).subs)
            %                             locs=locate_variables(s(1).subs,this.varnames,true);
            %                             if isnan(locs)
            %                                 % search within the date vector
            %                                 locs=locate_variables(s(1).subs,this.data(2:end,1),true);
            %                                 if isnan(locs)
            %                                     error([s(1).subs,' was not recognized either as a variable name or as a date'])
            %                                 end
            %                                 tmp=double(this);
            %
            %                             end
            %                         else
            %                             % locate among dates
            %                         end
        else%if  strcmp(s(1).type,'{}')
            access=true;
            if isnumeric(s(1).subs{1})
                access=false;
                if ~isscalar(s(1).subs{1})
                    error([mfilename,':: for leads or lags, subsref must be a scalar'])
                end
            elseif ischar(s(1).subs{1})||iscellstr(s(1).subs{1})
                % is it dates or is it varnames?
                if ischar(s(1).subs{1})
                    items=cellstr(s(1).subs{1});
                else
                    items=s(1).subs{1};
                end
                locs=locate_variables(items,this.varnames,true);
                by_cols=all(~isnan(locs));
                if ~by_cols
                    locs=locate_variables(items,this.data(2:end,1),true);
                    if any(isnan(locs))
                        error([mfilename,':: variables not found in database or dates not found in date vector'])
                    end
                end
            else
                error([mfilename,':: failed to understand this kind of subsref. Please contact junior.maih@gmail.com'])
            end
            datta=double(this);
            if access
                if by_cols
                    this=rise_time_series(this.start,datta(:,locs,:),this.varnames(locs));
                else
                    this=datta(locs,:,:);
                end
            else
                lead_or_lag=abs(s(1).subs{1});
                if s(1).subs{1}<0
                    datta=datta(1:end-lead_or_lag,:,:);
                    dates=serial2date(this.date_number(1+lead_or_lag));
                elseif s(1).subs{1}>0
                    datta=datta(1+lead_or_lag:end,:,:);
                    dates=serial2date(this.date_number(1));
                end
                this=rise_time_series(dates,datta,this.varnames);
            end
            s = s(2:end); % Remove processed s from list.
        end
    else
        this=subsref(this,s);
        s=[];
    end
end
end
