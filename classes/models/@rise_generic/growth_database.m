function db=growth_database(obj,endhist_date,end_sample,growth_type)

if isempty(obj)
    db=struct();
    return
end

if ~ismember(growth_type,{'zero','steady'})
    error('growth_type should be "zero" or "steady"')
end
zero_growth_=strcmp(growth_type,'zero');

% needs to be adapted for VARs
sv=get(obj,'state_vars');
maxlag=max(cell2mat(sv(2,:)))+1;
is_log_var=obj.endogenous.is_log_var;

ss=get(obj,'sstate');
g=get(obj,'growth');

endhist_date=date2serial(endhist_date);
starthist_date=endhist_date-maxlag+1;
last_date=date2serial(end_sample);
nperiods=numel(starthist_date:last_date);

db=struct();
dlin=zeros(nperiods,1);
dlog=ones(nperiods,1);
nv=obj.endogenous.number;
span=(1-maxlag:nperiods-maxlag)';
for iv=1:nv
    name=obj.endogenous.name{iv};
    if zero_growth_
        d=zero_growth();
    else
        d=steady_growth();
    end
    db.(name)=ts(starthist_date,d);
end

d=dlin;
for ix=1:sum(obj.exogenous.number)
    name=obj.exogenous.name{ix};
    db.(name)=ts(starthist_date,d);
end

    function d=zero_growth()
        if is_log_var(iv)
            d=dlog;
        else
            d=dlin;
        end
    end

    function d=steady_growth()
        ssv=ss.(name);
        gv=g.(name);        
        if is_log_var(iv)
            d=ssv*gv.^span;
        else
            d=ssv+gv*span;
        end
    end

end