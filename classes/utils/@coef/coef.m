classdef coef
    properties
        par_name
        func
        args
    end
    methods
        function obj=coef(varargin)
            % string
            % eqtn,var_pos,lag
            % eqtn,var_pos,lag,chain_name,state
            n=nargin;
            if n
                if isa(varargin{1},mfilename)
                    obj=varargin{1};
                else
                    if n==1
                        if ~ischar(varargin{1})
                            error('with one argument, coef expects a string or a coef')
                        end
                        obj.par_name=varargin{1};
                    else
                        obj.par_name=varargin;
                    end
                end
            end
        end
        function obj=plus(a,b)
            if ~(isa(a,'coef') && isa(b,'coef'))
                error('both arguments should be coef')
            end
            obj=coef();
            obj.args={a,b};
            obj.func='plus';
        end
        function obj=minus(a,b)
            if ~(isa(a,'coef') && isa(b,'coef'))
                error('both arguments should be coef')
            end
            obj=coef();
            obj.args={a,b};
            obj.func='minus';
        end
        function obj=mtimes(a,b)
            obj=coef();
            if isa(a,'coef') && isa(b,'coef')
                error('multiplication of coef x coef objects not allowed')
            end
            obj.args={a,b};
            obj.func='mtimes';
        end
        function obj=mrdivide(a,b)
            if ~isa(b,'double')
                error('the second argument must be a double')
            end
            obj=mtimes(a,1/b);
        end
        function [R,r]=linear_restrictions(obj,r,self,est_list)
            %--------------------------------
            [param_names,governing_chain,parameter_values,chain_names,...
                grand_chains_to_small,regimes,endo_names]=...
                main_object_environment();
            %--------------------------------
            nest=numel(est_list);
            nrest=numel(obj);
            if ~isequal(size(r),[nrest,1])
                error('RHS of restrictions (r) inconsitent with the number of restrictions')
            end
            R=sparse(nrest,nest);
            for iobj=1:nrest
                process_object(obj(iobj),iobj);
            end
            function idval=process_object(this,restr_id)
                idval=[];
                if isempty(this.args)
                    [id_,parval]=column_location(this.par_name);
                    if isempty(id_)
                        % elements to go to the rhs
                        if isnan(parval)
                            error([this.par_name,' is not estimated but found to be nan... this is not allowed'])
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
                    % separate the elements to go to the left from those
                    % that go to the right
                    rhs=imag(idval)==0;
                    idval_lhs=idval(~rhs);
                    idval_rhs=idval(rhs);
                    R(restr_id,real(idval_lhs))=imag(idval_lhs);
                    r(restr_id)=r(restr_id)-sum(idval_rhs);
                end
            end
            function [id,parval]=column_location(vv)
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
                        eqtn=vv{1};
                        var_pos=vv{2};
                        if ischar(var_pos)
                            var_pos=find(strcmp(var_pos,endo_names));
                            if isempty(var_pos)
                                error([var_pos,' is not recognized as a variable name'])
                            end
                        end
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
                    vv=parser.param_name_to_valid_param_name(vv);
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
    end
end