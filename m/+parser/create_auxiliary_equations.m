function [dic,modelBlock]=create_auxiliary_equations(dic,modelBlock)

[dic,modelBlock,LeadListing,new_auxvars,ilist_lead]=auxiliarize_exo_params(...
    dic,modelBlock);

[dic,modelBlock,LeadListing,new_auxvars,ilist_lead]=getRidOfExcessLeads(...
    dic,modelBlock,LeadListing,new_auxvars,ilist_lead);

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

function varargout=getRidOfExcessLags(varargin)
% function [dic,modelBlock,listing,new_auxvars,ilist]=getRidOfExcessLags(dic,...
%     modelBlock,listing,new_auxvars,ilist)

[varargout{1:nargout}]=getRidOfExcess(varargin{:},'-');

return

newMaxLeadOrLag=-1;

for irow=1:size(modelBlock,1)
    
    maxLag=modelBlock{irow,2};
    
    if maxLag>=-1
        
        continue
        
    end
    
    eqtn=modelBlock{irow,1};
    
    endog=dic.endogenous;
        
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
                
                newthing=sprintf('%s = %s{-1};',newguy,oldguy);
                
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

function varargout=getRidOfExcessLeads(varargin)
% function [dic,modelBlock,listing,new_auxvars,ilist]=getRidOfExcessLeads(...
%     dic,modelBlock,listing,new_auxvars,ilist)

[varargout{1:nargout}]=getRidOfExcess(varargin{:},'+');

return

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
                
                new_aux_eqtn=sprintf('%s{+1}=0;',old_guy);
                
            else
                
                newguy=makeNewJensenVariable(irow,start_index);
                
                new_aux_eqtn=sprintf('%s=%s{+1};',newguy,old_guy);
                
                % All new variables have a lead
                dic.endogenous=parser.update_variable_lead_lag(dic.endogenous,...
                    newguy,newMaxLeadOrLag,is_log_var);
                
                new_auxvars{ilist+2}=newguy;
                
                old_guy=newguy;
                
            end
            
            ilist=ilist+1;
            
            listing(ilist,:)={nan,new_aux_eqtn,'auxiliary equations'};
            
        end
        
    end

    function lagEquation()
        
        eqtn(2,:)=cellfun(@(x)if_then_else(~isempty(x),x-maxLead+1,x),...
            eqtn(2,:),'uniformOutput',false);
        
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


function [dic,modelBlock,listing,new_auxvars,ilist]=auxiliarize_exo_params(dic,modelBlock)

xo_param_names=[{dic.exogenous.name},{dic.parameters.name}];

ilist=0;

listing=cell(1000,3);

new_auxvars=cell(1,1000);

for irow=1:size(modelBlock,1)
    
    maxLag=modelBlock{irow,2};
    
    maxLead=modelBlock{irow,3};
    
    force_change=max(abs([maxLag,maxLead]))>1;
    
    [modelBlock{irow,1}]=auxiliarize_one_equation(...
        modelBlock{irow,1});
    
end

    function eqtn=auxiliarize_one_equation(eqtn)
        
        for icol0=1:size(eqtn,2)
            
            leadOrLag=eqtn{2,icol0};
            
            vname=eqtn{1,icol0};
            
            if ~any(strcmp(xo_param_names,vname))|| isempty(leadOrLag)
                
                continue
                
            end
            
            if ~force_change
                
                if leadOrLag==0
                    
                    continue
                    
                end
                
            end
            
%             if leadOrLag<0
%                 
%                 error('Leads or lags on parameters or exogenous not allowed')
%                 
%             end
            
            [new_name]=create_new_auxiliary_for_exo_or_param(vname,leadOrLag);
            
            eqtn{1,icol0}=new_name;
            
        end
        
        function  [new_name]=create_new_auxiliary_for_exo_or_param(...
                old_name,leadOrLag)
            
            endo_vars={dic.endogenous.name};
            
            new_name=parser.create_auxiliary_name(old_name,0,true);
            
            % if the name is not in the dictionary, create an auxiliary equation at the
            % same time
            
            if ~any(strcmp(endo_vars,new_name))
                
                ilist=ilist+1;
                
                new_aux_eqtn=sprintf('%s=%s;',new_name,old_name);
                
                listing(ilist,:)={nan,new_aux_eqtn,'auxiliary equations'};
                
                is_log_var=false;
                
                dic.endogenous=parser.update_variable_lead_lag(dic.endogenous,new_name,...
                    leadOrLag,is_log_var);
                
                new_auxvars{ilist}=new_name;
                
            end
            
        end

    end

end

function [dic,modelBlock,listing,new_auxvars,ilist]=getRidOfExcess(dic,...
    modelBlock,listing,new_auxvars,ilist,the_type)

switch the_type
    
    case '+'
        
        newMaxLeadOrLag=1;
        
        evalfunc=@le;
        
        update=@(x)x-1;
        
        theOne=1;
        
        relev_col=3;
        
    case '-'
        
        newMaxLeadOrLag=-1;
        
        evalfunc=@ge;
        
        update=@(x)x+1;
        
        theOne=-1;
        
        relev_col=2;

end

for irow=1:size(modelBlock,1)
    
    maxLag=modelBlock{irow,relev_col};
    
    if evalfunc(maxLag,newMaxLeadOrLag)
        
        continue
        
    end
    
    eqtn=modelBlock{irow,1};
    
    endog=dic.endogenous;
        
    for icol0=1:size(eqtn,2)
        
        theLag=eqtn{2,icol0};
        
        if isempty(theLag)||evalfunc(theLag,newMaxLeadOrLag)
            
            continue
            
        end
        
        endo_names={endog.name};
        
        vname=eqtn{1,icol0};
        
        newVname=parser.create_auxiliary_name(vname,update(theLag),true);
        
        eqtn{1,icol0}=newVname;
        
        eqtn{2,icol0}=newMaxLeadOrLag;
        
        pos=strcmp(vname,endo_names);
        
        is_log_var=endog(pos).is_log_var;
        
        endog(pos).max_lag=newMaxLeadOrLag;
        
        oldguy=vname;
        
        %----------------------------------------
        for item=1:abs(theLag)-1
            
            newguy=parser.create_auxiliary_name(vname,theOne*item,true);
            
            if ~any(strcmp(newguy,endo_names))
                
                endog=parser.update_variable_lead_lag(endog,newguy,theOne,...
                    is_log_var,vname);
                
                newthing=sprintf('%s = %s{%s1};',newguy,oldguy,the_type);
                
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
