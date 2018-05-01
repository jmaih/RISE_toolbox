function db=growth_database(obj,endhist_date,end_sample,growth_type)
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

if isempty(obj)

    db=cell(0,4);

    return

end

% needs to be adapted for VARs
sv=get(obj,'state_vars');

maxlag=max(cell2mat(sv(2,:)))+1;

endhist_date=date2serial(endhist_date);

starthist_date=endhist_date-maxlag+1;

last_date=date2serial(end_sample);

nperiods=numel(starthist_date:last_date);

span=(1-maxlag:nperiods-maxlag)';

endo_names=obj.endogenous.name;

nv=obj.endogenous.number;

exo_names=obj.exogenous.name;

nx=numel(exo_names);

db=struct();

regime=1;

if ischar(growth_type) 
    
    if ~(strncmp(growth_type,'zero',4)||...
            strncmp(growth_type,'steady',5))
        
        error('growth_type should be "zero" or "steady"')
        
    end
    
    growth_type(isspace(growth_type))=[];
    
    left_par=find(growth_type=='(');
    
    if ~isempty(left_par)
        
        right_par=find(growth_type==')');
        
        regime=str2double(growth_type(left_par+1:right_par-1));
        
        growth_type=growth_type(1:left_par-1);
        
    end
    
else
    
    if isa(growth_type,'ts')
        
        growth_type=pages2struct(growth_type);
        
    elseif ~isstruct(growth_type)
        
        error('growth_type must be "zero" or "steady" or a struct or a ts')
        
    end
    
    do_time_series();
    
    return
    
end

zero_growth_=strcmp(growth_type,'zero');

is_log_var=obj.endogenous.is_log_var;

ss=get(obj,'sstate');

g=get(obj,'growth');

d=[];

if zero_growth_
    
    span(:)=0;
    
end

for iv=1:nv
    
    name=endo_names{iv};
    
    do_growth();
            
    db.(name)=ts(starthist_date,d);
    
end

d=zeros(nperiods,1);

for ix=1:nx
    
    name=exo_names{ix};
    
    db.(name)=ts(starthist_date,d);
    
end

    function do_growth()
        
        ssv=ss.(name);
        
        gv=g.(name);
        
        if is_log_var(iv)
            
            d=ssv(regime)*gv(regime).^span;
            
        else
            
            d=ssv(regime)+gv(regime)*span;
            
        end
        
        d=full(d);
        
    end

    function do_time_series()
        
        growth_type=ts.collect(growth_type);
        
        varnames=growth_type.varnames;
        
        NumberOfPages=growth_type.NumberOfPages;
        
        NumberOfObservations=growth_type.NumberOfObservations;
        
        data=growth_type.data;
        
        allvars=[endo_names,exo_names];
        
        bad=find(~ismember(allvars,varnames));
        
            dn=growth_type.date_numbers;
            
        if any(bad)
            
            bang=find(date2serial(endhist_date)==dn);
            
            nbad=numel(bad);
            
            missing_names=allvars(bad);
            
            varnames=[varnames,missing_names];
            
            newdata=nan(NumberOfObservations,nbad,NumberOfPages);
            
            isexo=ismember(missing_names,exo_names);
            
            % set historical shocks to zero
            %-------------------------------
            newdata(1:bang,isexo)=0;
            
            data=cat(2,data,newdata);
            
            
        end
        
        is_start=find(starthist_date==dn);
        
        is_end=find(last_date==dn);
        
        data=data(is_start:is_end,:,:);
        
        db=ts(starthist_date,data,varnames);
        
        db=pages2struct(db);
    end

end