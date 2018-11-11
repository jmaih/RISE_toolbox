function [dic,modelBlock]=create_auxiliary_equations(dic,modelBlock)
% INTERNAL FUNCTION
%

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

end


function varargout=getRidOfExcessLeads(varargin)
% function [dic,modelBlock,listing,new_auxvars,ilist]=getRidOfExcessLeads(...
%     dic,modelBlock,listing,new_auxvars,ilist)

[varargout{1:nargout}]=getRidOfExcess(varargin{:},'+');

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
