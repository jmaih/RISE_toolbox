function [obj,retcode]=solve(obj,varargin)
% solver options are 0 or 'msre_klein' (for constant-parameter models)
%                    1 or 'msre_gensys' (for constant-parameter models)
%                    2 or 'msre_aim' (for constant-parameter models)
%                    3 or 'functional_iteration'
%                    4 or 'newton_kronecker'
%                    5 or 'newton_system'
%                    6 or 'newton_kronecker_iteration'
%                    7 or 'loose_com' (for constant-parameter optimal policy models)
% the algorithm automatically switches to the constant parameter solver
% 0 (Klein) whenever possible

Defaults=mergestructures(dsge_lc_solve(),msre_solve(),...
    struct('solver','msre_klein',...
    'check_stability',true,...
    'vectorized_code_threshold',150,...
    'accelerate_solver',true));

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    obj=Defaults;
    return
end

nobj=numel(obj);
if nobj>1
    retcode=nan(1,nobj);
    for iobj=1:nobj
        [obj(iobj),retcode(iobj)]=solve(obj(iobj),varargin{:});
    end
    return
end

params=[];
if ~isempty(varargin)
	tmp=find(strcmp('evaluate_params',varargin));
	if ~isempty(tmp)
		params=varargin{tmp+1};
		varargin(tmp+[0,1])=[];
	end
end
if ischar(params) && ~ismember(params,{'mode','mean'})
    error([mfilename,':: unrecognized string ''',params,''' for parameter input'])
end

n=numel(obj);
retcode=nan(1,n);
for jj=1:n
    if ~isempty(varargin)
        obj(jj).options=mysetfield(obj(jj).options,varargin{:});
    end
end
% varargin=rmfield(obj(1).options,...
%     setdiff(fieldnames(obj(1).options),fieldnames(Defaults)));
for jj=1:n
    [obj(jj),retcode(jj)]=solve_intern(obj(jj),params);
    if retcode(jj) && obj(jj).options.debug
        decipher_error(retcode(jj));
    end
end

    function [obj,retcode]=solve_intern(obj,params)
        retcode=0;
        
        % forward expansion order
        % if order==0, it means there is no ancitipation. It does not mean
        % that there no shocks. If order==0, this will later on help set
        % the conditional information to nan. But for the purpose of
        % expanding the model, the expansion order is at least 1.
        order=max(1,max(obj.options.irf_horizon,obj.options.order));
        forward_expansion_has_changed=~isequal(obj.current_expansion_order,order);
        solution_algorithm_has_changed=~isequal(obj.current_solution_algorithm,obj.options.solver);
        
        resolve_it=isempty(obj.T)||obj.estimation_under_way ||...
            (~isempty(params) && ~isequal(obj.current_solution_parameters,params))||... % ~isequal(obj.current_solution_parameters,params)||...
            forward_expansion_has_changed||...
            solution_algorithm_has_changed||...
			strcmp(obj.options.solve_initialization,'random');

        if forward_expansion_has_changed
            obj.current_expansion_order=order;
        end
        if solution_algorithm_has_changed
            obj.current_solution_algorithm=obj.options.solver;
        end
        
        
        if resolve_it
            [obj,retcode]=evaluate(obj,'evaluate_params',params);
            
            if ~retcode
                endo_nbr=obj.NumberOfEndogenous(2);
                h=obj.NumberOfRegimes;
                Aminus=obj.Aminus; A0=obj.A0; Aplus=obj.Aplus;
                B=obj.B; Q=obj.Q; hidden_ss_bgp=obj.steady_state_and_balanced_growth_path;
                
                T0=[]; % endo_nbr x endo_nbr x h initial guess for T
                
                % default solution for risk
                gsig=zeros(endo_nbr,h);
                
                if obj.is_hybrid_expectations_model
                    lambda=[obj.parameters(obj.hybrid_expectations_lambda_id).startval];
                    [Aplus,Aminus]=hybrid_expectations(Aplus,Aminus,lambda);
                end
                
                if obj.is_sticky_information_model
                    lambda=[obj.parameters(obj.sticky_information_lambda_id).startval];
                    % the hidden steady state is computed in evaluate
                    [Aplus,A0,Aminus,B,hidden_ss_bgp]=sticky_information(...
                        Aplus,A0,Aminus,B,hidden_ss_bgp,lambda,...
                        obj.forward_looking_ids);
                end
                
                if obj.is_optimal_policy_model
                    for ii=2:size(A0,3)
                        if ~isequal(A0(:,:,1),A0(:,:,ii))||~isequal(Aminus(:,:,1),Aminus(:,:,ii))||...
                                ~isequal(Aplus(:,:,1),Aplus(:,:,ii))||~isequal(B(:,:,1),B(:,:,ii))||...
                                ~isequal(obj.W(:,:,1),obj.W(:,:,ii))
                            error([mfilename,':: No switch allowed in the parameters of the loose commitment model'])
                        end
                    end
                    % do not take this steady state into consideration as
                    % it automatically returns 0. the user might have given
                    % a steady state file that returns different values.
                    % Plus it does not account for balanced growth, which
                    % is taken care of in the steady state computation.
                    [TT,RR,~,retcode,obj.options]=dsge_lc_solve(Aminus(:,:,1),A0(:,:,1),...
                        Aplus(:,:,1),B(:,:,1),obj.W(:,:,1),obj.planner_commitment(1),...
                        obj.planner_discount(1),order,obj.reordering_index,T0,obj.options);
                    % change the name of the solver right here, right now!
                    obj.options.solver='loose_com';
                    if ~retcode
                        [TT,RR,hidden_ss_bgp]=...
                            recast_loose_commitment_solution_into_markov_switching(obj,TT,RR,hidden_ss_bgp);
                    end
                else
                        loc=strcmp('is_switching',obj.parameters_image(1,:));
                        switching_params=obj.parameters_image{2,loc};
                        loc=strcmp('startval',obj.parameters_image(1,:));
                        sparam=obj.parameters_image{2,loc}(switching_params,:);
%                     [TT,RR,gsig,retcode,obj.options] = msre_solve_accelerated(Aminus,A0,Aplus,B,...
%                         obj.C,sparam,obj.is_unique_steady_state,...
%                         Q,order,T0,obj.options);
                    
                    [TT,RR,gsig,retcode,obj.options] = msre_solve(Aminus,A0,Aplus,B,...
                        obj.C,sparam,obj.is_unique_steady_state,...
                        Q,order,T0,obj.options);
                    
                    if obj.is_sticky_information_model
                        TT=TT(obj.reordering_index,obj.reordering_index,:);
                        RR=RR(obj.reordering_index,:,:,:);
                        hidden_ss_bgp=hidden_ss_bgp(obj.reordering_index,:,:);
                        gsig=gsig(obj.reordering_index,:); % this is unnecessary since all are zero (no parameter switch), but for completeness
                    end
                end
                
                if ~retcode
                    oldT=obj.T;
                    obj.T=TT;
                    if obj.options.check_stability && obj.NumberOfRegimes>1
                        if (ischar(obj.options.solver)&& ...
                                ~ismember(obj.options.solver,{'msre_gensys','msre_aim'}))||...
                                (isnumeric(obj.options.solver)&& ...
                                ~ismember(obj.options.solver,[1,2]))
                            if ~obj.is_stable_system
                                obj.T=oldT;
                                retcode=25; % system unstable
                            end
                        end
                    end
                    if ~retcode
                        % This is where I should apply the shock restrictions before passing everything
                        % to obj.R... The format of R is endo_nbr x exo_nbr x order x regime
                        if order>1 && ~isempty(obj.options.shock_properties)
                            nshocks=obj.NumberOfExogenous;
                            for ii=1:nshocks
                                shock=obj.options.shock_properties(ii).name;
                                shock_id=strcmp(shock,{obj.varexo.name});
                                horizon=obj.options.shock_properties(ii).horizon;
                                RR(:,shock_id,horizon+1:end,:)=0;
                            end
                        end
                        obj.R=RR;
                        obj.steady_state_and_balanced_growth_path=hidden_ss_bgp;
                        tmp=mat2cell(hidden_ss_bgp,ones(2*endo_nbr,1),h);
                        tmp0=tmp(1:endo_nbr,:);
                        % this thing below works well... even after checking many
                        % times :-)
                        [obj.varendo(:).det_steady_state]=(tmp0{:});
                        tmp0=tmp(endo_nbr+1:end,:);
                        [obj.varendo(:).balanced_growth]=(tmp0{:});
                        
                        tmp=mat2cell(gsig,ones(endo_nbr,1),h);
                        [obj.varendo(:).risk]=(tmp{:});
                        % store the following information so that we don't have to
                        % resolve the model if we don't have to. I probably should put
                        % an extra option called current_solver so that the model is
                        % resolved if the solver changes.
                        obj.current_solution_parameters=params;
                    end
                end
                %	% Sparsing does not work for arrays of more than 2 dimensions
                %	SparseList={'R','T'};
                %	for kk=1:numel(SparseList)
                %		obj.(SparseList{kk})=sparse(obj.(SparseList{kk}));
                %	end
                
            end
        end
    end
end

function [Aplus_he,Aminus_he]=hybrid_expectations(Aplus,Aminus,lambda)
if isscalar(lambda)
    Aminus_he=(1-lambda)*Aplus+Aminus;
    Aplus_he=lambda*Aplus;
else
    [eq_nbr,endo_nbr,regime_nbr]=size(Aplus);
    if numel(lambda)~=regime_nbr
        error([mfilename,':: lambda expected to have ',int2str(regime_nbr),' elements'])
    end
    Aminus_he=nan(eq_nbr,endo_nbr,regime_nbr);
    Aplus_he=nan(eq_nbr,endo_nbr,regime_nbr);
    for ii=1:regime_nbr
        Aminus_he(:,:,ii)=(1-lambda(ii))*Aplus(:,:,ii)+Aminus(:,:,ii);
        Aplus_he(:,:,ii)=lambda(ii)*Aplus(:,:,ii);
    end
end
end

function [Aplus_si,A0_si,Aminus_si,B_si,SS_and_BGP_si]=...
    sticky_information(Aplus,A0,Aminus,B,SS_and_BGP,lambda,truly_forward_id)
[old_endo_nbr,exo_nbr,regime_nbr]=size(B);
endo_nbr=old_endo_nbr+numel(truly_forward_id);
% initialize output matrices
Aplus_si=zeros(endo_nbr,endo_nbr,regime_nbr);
A0_si=repmat(eye(endo_nbr),[1,1,regime_nbr]);
Aminus_si=zeros(endo_nbr,endo_nbr,regime_nbr);
B_si=zeros(endo_nbr,exo_nbr,regime_nbr);
SS_and_BGP_si=zeros(2*endo_nbr,regime_nbr);
% define the steady_state and balance growth paths locations
ss=1:old_endo_nbr; bgp=old_endo_nbr+1:2*old_endo_nbr;
fw=numel(truly_forward_id);
orig_endo_ss_and_bgp=[ss,bgp+fw];
sticky_ss_and_bgp=[ss+(1:fw),endo_nbr+old_endo_nbr+(1:fw)];

% the variables are arranged as : [Y',SI']' with the sticky information
% variables below, Y_{t+1}=lamb*Y_{t+1}+(1-lamb)*SI_{t-1} and
% SI_t=E_t(Y_{t+1})
B_si(1:old_endo_nbr,:,:)=B;
SS_and_BGP_si(orig_endo_ss_and_bgp,:)=SS_and_BGP;
A0_si(1:old_endo_nbr,1:old_endo_nbr,:)=A0;
Aminus_si(1:old_endo_nbr,1:old_endo_nbr,:)=Aminus;
% add the steady states and BGP for the sticky information variables
SS_and_BGP_si(sticky_ss_and_bgp,:)=...
    SS_and_BGP([ss(truly_forward_id),bgp(truly_forward_id)],:);

if isscalar(lambda)
    Aminus_si(1:old_endo_nbr,old_endo_nbr+1:end,:)=(1-lambda)*Aplus(:,truly_forward_id,:);
    Aplus_si(1:old_endo_nbr,1:old_endo_nbr,:)=lambda*Aplus;
else
    if numel(lambda)~=regime_nbr
        error([mfilename,':: lambda expected to have ',int2str(regime_nbr),' elements'])
    end
    for ii=1:regime_nbr
        Aminus_si(1:old_endo_nbr,old_endo_nbr+1:end,ii)=(1-lambda(ii))*Aplus(:,truly_forward_id,ii);
        Aplus_si(1:old_endo_nbr,1:old_endo_nbr,ii)=lambda(ii)*Aplus(:,:,ii);
    end
end
Aplus_si(old_endo_nbr+1:end,truly_forward_id,:)=...
    repmat(-eye(numel(truly_forward_id)),[1,1,regime_nbr]);
end