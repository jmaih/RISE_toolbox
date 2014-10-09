function this=subsasgn(this,s,b)% subsasgn
if strcmp(s(1).type,'.') || length(s)>1
    this=builtin('subsasgn',this,s,b);
else
    vcol=1;
    rowsid=1;
    page_id=1;
    for ii=1:length(s.subs)
        current_subs=s.subs{ii};
        % load the serial number if it is date
        serial=[];
        if ischar(current_subs)||isnumeric
            % check whether it is a date
            serial=date2serial(current_subs,true); % true for not throwing an error inside date2serial
        end
        if ~isempty(serial)
            nserial=numel(serial);
            rowsid=nan(nserial,1);%
            for irow=1:nserial
                tmp=find(this.date_number==serial(irow));
                if ~isempty(tmp)
                    rowsid(irow)=tmp;
                end
            end
            if any(isnan(rowsid))
                error('the dates above could not be located')
            end
        elseif ischar(current_subs)
            loc=find(strcmp(current_subs,this.varnames));
            if ~isempty(loc)
                vcol=loc;
            else
                error([mfilename,':: could not identify string',current_subs])
            end
        elseif isnumeric(current_subs) && ismember(current_subs,[1,this.NumberOfPages])
            page_id=current_subs;
        else
            error([mfilename,':: I do not understand the assignment'])
        end
    end
    datta=double(this);
    if ~isscalar(b) && (numel(b)~=numel(rowsid))
        error([mfilename,':: number of rhs elements does not match lhs'])
    end
    datta(rowsid,vcol,page_id)=b;
    this=rise_time_series(this.start,datta,this.varnames);
end
end
