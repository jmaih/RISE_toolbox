function [dictionary,PlannerObjective_block,Model_block]=planner_objective(...
dictionary,listing,Model_block)
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

dictionary.is_model_with_planner_objective=~isempty(listing);

PlannerObjective_block=cell(0,1);

dictionary.planner_system=[];

if ~dictionary.is_model_with_planner_objective
    
    return

end

% there can't be multiple planner_objective blocks
% pull everything on one line
fichier_name_=listing{1,3};

line_number_=listing{1,1};

OneLiner='';

for ii=1:size(listing,1)
    
    OneLiner=[OneLiner,listing{ii,2}]; %#ok<*AGROW>
    
    if ~ismember(listing{ii,1},line_number_)
        
        line_number_=[line_number_,listing{ii,1}];
    
    end
    
end

OneLiner(isspace(OneLiner))=[];

Commitment='commitment=1;';

Discount='discount=.99;';

Objective='';

[startd,endd]=regexp(OneLiner,'discount=\w*\.?\w+','start','end');

if ~isempty(startd)
    
    Discount=[OneLiner(startd:endd),';'];
    
    OneLiner(startd:endd)=[];

end

[startc,endc]=regexp(OneLiner,'commitment=\w*\.?\w+','start','end');

if ~isempty(startc)
    
    Commitment=[OneLiner(startc:endc),';'];
    
    OneLiner(startc:endc)=[];
    
    % Time to decide whether to add a markov chain controlling policy
    % behavior or not. This depends on whether Commitment is loose
    equality=strfind(Commitment,'=');
    
    true_commitment=Commitment(equality+1:end-1);
    
    if ~ismember(true_commitment,{'0','1'})
        
        dictionary.markov_chains(end+1)=parser.initialize_markov_chain(...
            parser.loose_commit(),2,'is_endogenous',false);
        
        % Adding a new markov chain through loose commitment. In
        % re-ordering the chains, we need to update the info in the
        % parameters
        [~,tmp]=sort({dictionary.markov_chains.name});
        
        for iparam=1:numel(dictionary.parameters)
            
            dictionary.parameters(iparam).governing_chain=...
                find(dictionary.parameters(iparam).governing_chain==tmp);
        
        end
        
        dictionary.markov_chains=dictionary.markov_chains(tmp);
    
    end
    
end

if ~isempty(OneLiner)
    
    right_curly_brace=strfind(OneLiner,'}');
    
    if ~isempty(right_curly_brace)
        
        OneLiner=OneLiner(right_curly_brace+1:end);
    
    end
    
    Objective=OneLiner;

end

if strcmp(Objective(1),';')
    
    error([mfilename,':: No objective function provided. Check ',fichier_name_,' at line(s) ',mat2str(line_number_)])

end

tmp={line_number_,Objective,fichier_name_; % policy objective
    line_number_,Commitment,fichier_name_; % commitment
    line_number_,Discount,fichier_name_};  % discount

[PlannerObjective_block,dictionary]=parser.capture_equations(dictionary,tmp,'planner_objective');

PlannerObjective_block(1,:)=remove_leads_and_lags(PlannerObjective_block(1,:));

dictionary.planner_system.model={Objective;Commitment;Discount};

    function in=remove_leads_and_lags(in)
        
        time_vars=cell2mat(in(2:3));
        
        if isempty(time_vars)
            
            error([mfilename,':: Planner objective must include variables. Check ',fichier_name_,' at line(s) ',mat2str(line_number_)])
            
        end
        
        if any(time_vars~=0)
            
            ilist=0;
            
            myListing=cell(50,3);
            
            new_auxvars=cell(1,50);
            
            eqtn=in{1};
            
            endog=dictionary.endogenous;
            
            for iii=1:size(in{1},2)
                
                lagOrLead=eqtn{2,iii};
                
                if isempty(lagOrLead)||(lagOrLead==0)
                    
                    continue
                    
                end
                
                endo_names={endog.name};
                
                vname=eqtn{1,iii};
                
                newVname=parser.create_auxiliary_name(vname,lagOrLead,true);
                
                eqtn{1,iii}=newVname;
                
                eqtn{2,iii}=0;
                
                pos=strcmp(vname,endo_names);
                
                is_log_var=endog(pos).is_log_var;
                
                endog(pos).max_lag=min(lagOrLead,endog(pos).max_lag);
                
                endog(pos).max_lead=max(lagOrLead,endog(pos).max_lead);
                
                if ~any(strcmp(newVname,endo_names))
                    
                    endog=parser.update_variable_lead_lag(endog,newVname,...
                        0,is_log_var,vname); % 0 is the new lag and lead for the new variable
                    
                    addstr='';
                    
                    if lagOrLead>0
                        
                        addstr='+';
                        
                    end
                    
                    the_type=sprintf('%s%0.0f',addstr,lagOrLead);
                    
                    newthing=sprintf('%s = %s{%s};',newVname,vname,the_type);
                    
                    ilist=ilist+1;
                    
                    myListing(ilist,:)={nan,newthing,'auxiliary equations'};
                    
                    new_auxvars{ilist}=newVname;
                    
                end
                
            end
            
            myListing=myListing(1:ilist,:);
            
            % push the modified equation
            %---------------------------
            in{1}=eqtn;
            
            % reset max leads and lags... or rather keep them as is to
            % remember what they were
            %-----------------------------------------------------------
            % in{2}=0; in{3}=0;
            
            dictionary.endogenous=endog;
            
            % add the equations to the model block
            [aux,dictionary]=parser.capture_equations(dictionary,myListing,'model');
            
            Model_block=[Model_block;aux];
                        
            % sweep once more to remove leads and lags that were pushed
            % into auxiliary equations...            
            [dictionary,Model_block]=parser.create_auxiliary_equations(dictionary,Model_block);
            
        end
        
    end

end

