function [R,r]=linear_restrictions(obj,r,self,est_list)
% This function construct the linear restrictions system R*a=r
% given the coef vector "obj", the "r" vector which may be
% modified, "self", which is a rise object (dsge, rfvar, svar,
% etc.) and "est_list", which is the list of the parameters to
% estimate.

% load the main object fields
%----------------------------
[param_names,governing_chain,parameter_values,chain_names,...
    grand_chains_to_small,regimes,endo_names]=...
    main_object_environment();
est_list=parser.param_texname_to_param_name(est_list);
nest=numel(est_list);
nrest=numel(obj);
if ~isequal(size(r),[nrest,1])
    error('RHS of restrictions (r) inconsitent with the number of restrictions')
end
R=sparse(nrest,nest);
for iobj=1:nrest
    process_object(obj(iobj),iobj);
end

% Now, locate and re-scale the dirichlet columns. This has to be done
% before any further processing of R
%--------------------------------------------------------------------------
correct_dirichlet()

    function correct_dirichlet()
        % Because the dirichlet will be unweighted i.e. xj-->xj*s. In order to
        % remain consistent with linear restrictions, the coefficients on the
        % dirichlet have to be divided by s. For that we use the location
        % information
        for id=1:numel(self.estim_dirichlet)
            di=self.estim_dirichlet(id);
            locs=di.location;
            h=numel(locs)+1;
            s=utils.distrib.dirichlet_sum_weights(h);
            R(:,locs)=R(:,locs)/s;
        end
    end

    function idval=process_object(this,restr_id)
        idval=[];
        if isa(this,'double')
            idval=this;
        elseif isempty(this.args)
            [id_,parval,pname]=column_location(this.par_name);
            if isempty(id_)
                % elements to go to the rhs
                if isnan(parval)
                    error([pname,' is not estimated but found to be nan... this is not allowed'])
                end
                idval=parval;
            else
                % elements to stay on the lhs
                idval=id_+1i;
            end
        else
            % go through the tree
            args_=this.args;
            if isa(args_{2},'double')
                args_=args_(end:-1:1);
            end
            switch this.func
                case 'plus'
                    idval1=process_object(args_{1});
                    idval2=process_object(args_{2});
                    idval=[idval1,idval2];
                case 'minus'
                    idval1=process_object(args_{1});
                    idval2=process_object(args_{2});
                    idval=[idval1,real(idval2)-1i*imag(idval2)];
                case 'mtimes'
                    idval=process_object(args_{2});
                    idval=real(idval)+args_{1}*imag(idval)*1i;
                otherwise
            end
        end
        if nargin>1 && ~isempty(idval)
            % separate the elements to go to the left from those that go to
            % the right
            rhs=imag(idval)==0;
            idval_lhs=idval(~rhs);
            idval_rhs=idval(rhs);
            R(restr_id,real(idval_lhs))=imag(idval_lhs);
            r(restr_id)=r(restr_id)-sum(idval_rhs);
        end
    end
    function [id,parval,vv]=column_location(vv)
        reformat_parameter_name();
        id=find(strcmp(vv,est_list));
        parval=nan;
        if isempty(id)
            % search the non-estimated list
            [vv_pos,vv_regime_states]=...
                generic_tools.parameter_position_and_regimes(vv,...
                param_names,governing_chain,chain_names,...
                grand_chains_to_small,regimes);
            parval=parameter_values(vv_pos,vv_regime_states{1}(1));
        end
        function reformat_parameter_name()
            if ~ischar(vv)% %eqtn,var_pos,lag,chain_name,state
                eqtn=find_variable_position(vv{1});
                var_pos=find_variable_position(vv{2});
                lag=vv{3};
                pname=sprintf('a%0.0f_%0.0f_%0.0f',...
                    lag,eqtn,var_pos);
                if numel(vv)>3
                    chain_name=vv{4};
                    state=vv{5};
                    pname=sprintf('%s(%s,%0.0f)',pname,chain_name,state);
                end
                vv=pname;
            end
            vv=parser.param_texname_to_param_name(vv);
            function loc=find_variable_position(str)
                loc=str;
                if ischar(str)
                    loc=find(strcmp(str,endo_names));
                    if isempty(str)
                        error([str,' is not recognized as a variable name'])
                    end
                end
                if loc>numel(endo_names)
                    disp(vv)
                    error(['One of the indexes in the name above ',...
                        'exceeds the number of endogenous variables'])
                end
            end
        end
    end
    function [param_names,governing_chain,parameter_values,...
            chain_names,grand_chains_to_small,regimes,endo_names]=...
            main_object_environment()
        parameter_values=self.parameter_values;
        param_names=self.parameters.name;
        grand_chains_to_small=self.markov_chains.grand_chains_to_small;
        par_nbr=sum(self.parameters.number);
        regimes=cell2mat(self.markov_chains.small_markov_chain_info.regimes(2:end,2:end));
        reg_nbr=size(regimes,1);
        if isempty(self.parameter_values)
            self.parameter_values=nan(par_nbr,reg_nbr);
        end
        chain_names=self.markov_chains.small_markov_chain_info.chain_names;
        governing_chain=self.parameters.governing_chain;
        endo_names=self.endogenous.name;
    end
end
