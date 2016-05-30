function [dic,modelBlock]=create_auxiliary_equations(dic,modelBlock)

[dic,modelBlock,LeadListing,new_auxvars,ilist_lead]=getRidOfExcessLeads(dic,modelBlock);

[dic,modelBlock,BigListing,new_auxvars,ilist_all]=getRidOfExcessLags(dic,...
    modelBlock,LeadListing,new_auxvars,ilist_lead); %#ok<ASGLU>

% add the equations to the model block
[aux,dic]=parser.capture_equations(dic,BigListing,'model');

modelBlock=[modelBlock;aux];

if isfield(dic,'auxiliary_equations')
    
    dic.auxiliary_equations=[dic.auxiliary_equations;aux];
    
else
    
    dic.auxiliary_equations=aux;
    
end

dic.auxiliary_variables.model=[dic.auxiliary_variables.model,...
    new_auxvars];

end

function [dic,modelBlock,listing,new_auxvars,ilist]=getRidOfExcessLags(dic,...
    modelBlock,listing,new_auxvars,ilist)

newMaxLeadOrLag=-1;

for irow=1:size(modelBlock,1)
    
    maxLag=modelBlock{irow,2};
    
    if maxLag>=-1
        
        continue
        
    end
    
    eqtn=modelBlock{irow,1};
    
    endog=dic.endogenous;
    
    str='-';
    
    for icol0=1:size(eqtn,2)
        
        theLag=eqtn{2,icol0};
        
        if isempty(theLag)||theLag>=-1
            
            continue
            
        end
        
        endo_names={endog.name};
        
        vname=eqtn{1,icol0};
        
        newVname=parser.create_auxiliary_name(vname,theLag+1,true);
        
        eqtn{1,icol0}=newVname;
        
        eqtn{2,icol0}=newMaxLeadOrLag;
        
        pos=strcmp(vname,endo_names);
        
        is_log_var=endog(pos).is_log_var;
        
        endog(pos).max_lag=newMaxLeadOrLag;
        
        oldguy=vname;
        
        %----------------------------------------
        for item=1:abs(theLag)-1
            
            newguy=parser.create_auxiliary_name(vname,-item,true);
            
            if ~any(strcmp(newguy,endo_names))
                
                endog=parser.update_variable_lead_lag(endog,newguy,-1,...
                    is_log_var,vname);
                
                newthing=sprintf('%s = %s{%s1};',newguy,oldguy,str);
                
                ilist=ilist+1;
                
                listing(ilist,:)={nan,newthing,'auxiliary equations'};
                
                endo_names={endog.name};
                
                new_auxvars{ilist}=newguy;
            end
            
            oldguy=newguy;
            
        end
        %----------------------------------------
        
    end
    
    dic.endogenous=endog;
    
    modelBlock{irow,1}=eqtn;
      
end

listing=listing(1:ilist,:);

new_auxvars=new_auxvars(1:ilist);

end

function [dic,modelBlock,listing,new_auxvars,ilist]=getRidOfExcessLeads(dic,modelBlock)

% Do not count these as auxiliary

ilist=0;

listing=cell(1000,3);

new_auxvars=cell(1,1000);

newMaxLeadOrLag=1;

is_log_var=false;

for irow=1:size(modelBlock,1)
    
    maxLead=modelBlock{irow,3};
    
    if maxLead<=1
        
        continue
        
    end
    
    eqtn=modelBlock{irow,1};
    
    start_index=0;
    
    auxVar=makeNewJensenVariable(irow,start_index);
    
    lagEquation()
    
    modelBlock{irow,1}=eqtn;
    
    createAuxiliaryLeads()
    
end

    function createAuxiliaryLeads()
        
        old_guy=auxVar;
        
        for ilead=2:maxLead
            
            start_index=start_index+1;
            
            if ilead==maxLead
                
                newthing=sprintf('%s{+1}=0;',old_guy);
                
            else
                
                newguy=makeNewJensenVariable(irow,start_index);
                
                newthing=sprintf('%s=%s{+1};',newguy,old_guy);
                
                % All new variables have a lead
                dic.endogenous=parser.update_variable_lead_lag(dic.endogenous,...
                    newguy,newMaxLeadOrLag,is_log_var);
                
                new_auxvars{ilist+2}=newguy;
                
                old_guy=newguy;
                
            end
            
            ilist=ilist+1;
            
            listing(ilist,:)={nan,newthing,'auxiliary equations'};
            
        end
        
    end

    function lagEquation()
        
        eqtn(2,:)=cellfun(@(x)if_then_else(~isempty(x),x-3,x),eqtn(2,:),...
            'uniformOutput',false);
        
        eqtn=[
            {'-',auxVar,'+'
            [],0,[]},eqtn];
        
        leadLags=cell2mat(eqtn(2,:));
        
        modelBlock{irow,3}=max(leadLags);
        
        leadLags(leadLags>0)=[];
        
        modelBlock{irow,2}=min(leadLags);
        
        % add variable to the dictionary: lead is 1 since the variable will be
        % lead shortly
        dic.endogenous=parser.update_variable_lead_lag(dic.endogenous,...
            auxVar,newMaxLeadOrLag,is_log_var);
        
        new_auxvars{ilist+1}=auxVar;
        
        % update leads and lags of all the variables in the equation
        %-----------------------------------------------------------
        for icol=1:size(eqtn,2)
            
            theLagOrLead=eqtn{2,icol};
            
            if isempty(theLagOrLead)
                
                continue
                
            end
            
            dic.endogenous=parser.update_variable_lead_lag(dic.endogenous,...
                eqtn{1,icol},theLagOrLead,is_log_var);
            
        end
        
    end

    function av=makeNewJensenVariable(irow,index)
        
        av=sprintf('EQTN_%0.0f_%0.0f',irow,index);
        
    end

end
