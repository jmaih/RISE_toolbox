function endo_struct=update_variable_lead_lag(endo_struct,vname,lead_lag,...
    is_log_var,original_name)

if nargin<5
    
    original_name=vname;
    
end

endo_names={endo_struct.name};

vpos=find(strcmp(vname,endo_names),1);

if ~isempty(vpos)
    
    new_var=endo_struct(vpos);
    
else
    
    vpos=numel(endo_names)+1;
    
    new_var=parser.listing('name',vname,'current_name',original_name,...
        'is_auxiliary',true,'is_log_var',is_log_var);
    
end

new_var.max_lead=max(new_var.max_lead,lead_lag);

new_var.max_lag=min(new_var.max_lag,lead_lag);

endo_struct(vpos)=new_var;

end