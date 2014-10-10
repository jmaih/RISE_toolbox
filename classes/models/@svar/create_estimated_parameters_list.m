function estim_names=create_estimated_parameters_list(obj)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 


if isempty(obj)
    obj=struct();
    return
end

% create baseline from the template for the estimation
%-----------------------------------------------------
[plist,fixed_struct]=vartools.create_baseline_parameters(obj.estim_param_template,obj.parameters.name);

Regimes=obj.markov_chains.regimes;
nstates=nan(1,obj.markov_chains.chains_number);
Regimes_mat=cell2mat(Regimes(2:end,2:end));
for istate=1:obj.markov_chains.chains_number
    nstates(istate)=max(Regimes_mat(:,istate));
end
chain_names=Regimes(1,2:end);
npar=sum(obj.parameters.number);

estim_names=generate_estim_names();

    function estim_names=generate_estim_names()
        incmnt=300;
        nrows2=incmnt;
        iter2=0;
        estim_names=cell(1,nrows2);
        for ip_=1:npar
            pname=obj.parameters.name{ip_};
            % transition probabilities are automatically added
            %-------------------------------------------------
            if ismember(pname,plist)||parser.is_transition_probability(pname)
                gov_chain=obj.parameters.governing_chain(ip_);
                if nstates(gov_chain)==1
                    store_name(pname);
                else
                    for istate_=1:nstates(gov_chain)
                        v=sprintf('%s_%s_%0.0f',pname,chain_names{gov_chain},istate_);
                        store_name(v);
                    end
                end
            else
                % set the parameter in all regimes the value in the template
                %----------------------------------------------------------
                obj.parameter_values(ip_,:)=fixed_struct.(pname);
            end
        end
        estim_names=estim_names(1:iter2);
        function store_name(v)
            iter2=iter2+1;
            if nrows2<iter2
                estim_names=[estim_names,cell(1,incmnt)];
                nrows2=nrows2+incmnt;
            end
            estim_names{iter2}=v;
        end
    end
end

% % if isa(obj,'stochvol')
% %     plist=setdiff(obj.parameters.name,fieldnames(fixed_struct));
% %     % taking just a setdiff with respect to the fixed is not enough
% %     % because the random walk assumption will also fix some of the
% %     % parameters
% %     % add the parameters of the time-varying parameter model
% %     %-------------------------------------------------------
% %     if obj.time_varying_parameters && obj.random_walk_parameters
% %         % fix b=0 and rho=1
% %         further_fix=regexp(plist,'(?<!\w+)a\d+_\d+_\d+(?!\w+)','match');
% %         further_fix=[further_fix{:}];
% %         fixed_struct=fix_parameters(fixed_struct,further_fix,0);
% %         further_fix=regexp(plist,'(?<!\w+)rho_a\d+_\d+_\d+(?!\w+)','match');
% %         further_fix=[further_fix{:}];
% %         fixed_struct=fix_parameters(fixed_struct,further_fix,1);
% %     end
% % % %     % add the parameters of the stochastic volatility restrictions
% % % %     %-------------------------------------------------------------
% % % %     further_fix=regexp(plist,'(?<!\w+)sig_\d+_\d+(?!\w+)','match');
% % % %     further_fix=[further_fix{:}];
% % % %     fixed_struct=fix_parameters(fixed_struct,further_fix,0);
% % end

%     function xx=fix_parameters(xx,fix_list,d)
%         for ilist=1:numel(fix_list)
%             xx.(fix_list{ilist})=d;
%         end
%         plist=setdiff(plist,fix_list);
%     end
