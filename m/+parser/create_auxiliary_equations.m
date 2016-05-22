function [dic,modelBlock]=create_auxiliary_equations(dic,modelBlock)

has_aux=isfield(dic,'auxiliary_equations') && ~isempty(dic.auxiliary_equations);

if has_aux
    
    old_aux=dic.auxiliary_equations(:,1);
    
    for iaux=1:numel(old_aux)
    
        old_aux{iaux}=old_aux{iaux}{1};

    end
    
else
    
    old_aux=cell(0,1);

end

listing=cell(100,3);

ilist=0;

endog=dic.endogenous;

endo_nbr=numel(endog);

new_auxvars=cell(1,100);

for ivar=1:endo_nbr
    
    this=endog(ivar);
    
    is_log_var=this.is_log_var;
    
    excess_lags=abs(this.max_lag);
    
    excess_leads=abs(this.max_lead);
    
    if max(excess_leads,excess_lags)>1
        
        vname=this.name;
        
        current_name=this.current_name;
        
        do_one_variable(excess_lags,-1);
    
        do_one_variable(excess_leads,1);

    end
    
end

dic.endogenous=endog;

listing=listing(1:ilist,:);

if isfield(dic,'latent_equations')
    % add the main auxiliary first
    %------------------------------
    listing=[dic.latent_equations;listing];
    % remove the field immediately
    %------------------------------
    dic=rmfield(dic,'latent_equations');

end

[aux,dic]=parser.capture_equations(dic,listing,'model');

modelBlock=[modelBlock;aux];

if ~isfield(dic,'auxiliary_equations')

    dic.auxiliary_equations=aux;

else
    
    dic.auxiliary_equations=[dic.auxiliary_equations;aux];

end

dic.auxiliary_variables.model=[dic.auxiliary_variables.model,...
    new_auxvars(1:ilist)];

    function do_one_variable(n,seign)
        
        oldguy=vname;
        
        if seign>0
        
            str='+';
        
        else
            
            str='-';
        
        end
        
        for item=1:n-1
            
            newguy=parser.create_auxiliary_name(vname,seign*item,true);
            
            endog=parser.update_variable_lead_lag(endog,newguy,seign,is_log_var,...
                current_name);
            
            newthing=sprintf('%s = %s{%s1};',newguy,oldguy,str);
            
            oldguy=newguy;
            
            if item==n-1
                
                if seign>0
                
                    endog(ivar).max_lead=1;
                
                else
                    
                    endog(ivar).max_lag=1;
            
                end
                
            end
            
            if ~any(strcmp(newguy,old_aux))
                % avoid duplicates
                ilist=ilist+1;
                
                listing(ilist,:)={nan,newthing,'auxiliary equations'};
            
                new_auxvars{ilist}=newguy;
            
            end
            
        end
        
    end

end