function [structural_matrices,retcode]=evaluate_all_derivatives(obj,structural_matrices,ssdata)
% evaluate all the derivatives up to the desired order: analytical
% derivatives can be evaluated sequentially, but for algorithmic
% derivatives and for numerical derivatives, it is more economical to do
% them at once, instead of calling the same functions several times

if nargin < 3
    
    ssdata = [];
    
end

solve_order=obj.options.solve_order;

h=obj.markov_chains.small_markov_chain_info.regimes_number;

nx=sum(obj.exogenous.number);

xss=zeros(nx,1);

params=obj.parameter_values;

def=obj.solution.definitions;

ss=cell2mat(obj.solution.ss);

is_has_data=~isempty(ssdata);

if ~is_has_data
    
    ssdata=ss;
    
end

if size(ssdata,2)~=h
    
    if size(ssdata,2)~=1
        
        error('wrong number of columns for ssdata')
        
    else
        
        ssdata=ssdata(:,ones(h,1));
        
    end
    
end

if size(ssdata,1)~=size(ss,1)
    
    error('wrong number of rows for ssdata')
    
end

[ys,nind]=forms_used_in_computation_of_derivatives();

if is_has_data
    % update residuals with the data
    for ireg=1:h
        
        s0=ireg;
        
        s1=ireg;
        
        structural_matrices.user_resids(:,ireg)=utils.code.evaluate_functions(...
            obj.routines.probs_times_dynamic,...
            ys(:,s0),xss,ss(:,s0),params(:,s0),def{s0},s0,s1);
        
    end
    
end

% evaluate first-order derivatives
%---------------------------------
symbolic_type=strcmp(obj.options.solve_derivatives_type,'symbolic');

automatic_type=false;

if ~symbolic_type
    
    automatic_type=strcmp(obj.options.solve_derivatives_type,'automatic');
    
    if ischar(obj.options.solve_automatic_differentiator)
        
        obj.options.solve_automatic_differentiator=...
            func2str(obj.options.solve_automatic_differentiator);
        
    end
    
    if ~automatic_type
        
        numeric_type=any(strcmp(obj.options.solve_derivatives_type,{'numeric','numerical'}));
        
        if ~numeric_type
            
            error(['solve_derivatives_type can only assume values ',...
                '"symbolic", "automatic" or "numeric"'])
            
        end
        
        if solve_order>1
            
            error('numerical derivatives not implemented for orders greater than 1')
            
        end
        % prepare the re-ordering of the endogenous columns
        reordering=obj.lead_lag_incidence.before_solve(obj.order_var,:);
        
        reordering=reordering(reordering>0);
        
    end
    
end

xxx=repmat('v',1,solve_order);
% evaluate higher-order derivatives
%----------------------------------
retcode=0;

G01=cell(1,solve_order);

for s1=1:h
    
    for s0=1:h
        
        if ~retcode
            % Note: G(s0,s1) =: ps0(s0,s1)*F(s0)
            if symbolic_type
                
                if s1==1 && s0==1
                    
                    max_order=numel(obj.routines.probs_times_dynamic_derivatives);
                    
                    if solve_order>max_order
                        
                        error(['Perturbation of order ',int2str(solve_order),...
                            ' requested but symbolic derivatives available ',...
                            'only up to order ',int2str(max_order),...
                            '. Compute higher-order symbolic derivatives or ',...
                            'switch to automatic derivatives'])
                        
                    end
                    
                end
                
                [G01{1:solve_order}]=utils.code.evaluate_functions(obj.routines.probs_times_dynamic_derivatives,...
                    ys(:,s0),xss,ss(:,s0),params(:,s0),def{s0},s0,s1);
                
            elseif automatic_type
                
                if s1==1 && s0==1
                    
                    max_order=sum(cell2mat(regexp(fieldnames(aplanar),'dx+')));
                    
                    if solve_order>max_order
                        
                        error(['Perturbation of order ',int2str(solve_order),...
                            ' requested but automatic derivatives available ',...
                            'only up to order ',int2str(max_order)])
                        
                    end
                    
                end
                
                engine=obj.options.solve_automatic_differentiator;
                
                G01=utils.code.evaluate_automatic_derivatives(...
                    obj.routines.symbolic.probs_times_dynamic,...
                    solve_order,engine,...
                    ys(:,s0),xss,ss(:,s0),params(:,s0),def{s0},s0,s1);
                
            else
                
                if s1==1 && s0==1
                    
                    max_order=1;
                    
                    if solve_order>max_order
                        
                        error(['Perturbation of order ',int2str(solve_order),...
                            ' requested but numerical derivatives available ',...
                            'only up to order ',int2str(max_order),...
                            '. Switch to symbolic or automatic derivatives'])
                        
                    end
                    
                end
                
                [G01{1:1}]=utils.code.evaluate_jacobian_numerically(obj.routines.probs_times_dynamic,...
                    ys(:,s0),xss,ss(:,s0),params(:,s0),def{s0},s0,s1);
                % The columns are to be put in the order_var order
                
                G01{1}(:,1:nind)=G01{1}(:,reordering);
                
            end
            
            if utils.error.valid(G01)
                % use the derivatives Gi to build dv, dvv, dvvv, ...
                %---------------------------------------------------
                zkz=1;
                
                log_deriv_coefs=structural_matrices.log_deriv_coefs.';
                
                for io=1:solve_order
                    
                    zkz=kron(zkz,log_deriv_coefs(s1,:));
                    structural_matrices.(['d',xxx(1:io)]){s0,s1}=...
                        bsxfun(@times,G01{io},zkz);%
                    
                end
                
            else
                
                retcode=2; % nans in jacobian
                
            end
            
        end
        
    end
    
end
% Compute planner information first
%----------------------------------
if obj.is_optimal_policy_model|| obj.is_optimal_simple_rule_model
    
    planner=struct(...
        'objective',{cell(1,h)},...
        'commitment',{cell(1,h)},...
        'discount',{cell(1,h)},...
        'weights',{cell(1,h)}...
        );
    
    orig_endo_nbr=obj.endogenous.number;
    
    for s0=1:h
        
        if ~retcode
            % approximation taken around the data and not the steady state
            %-------------------------------------------------------------
            lcd=utils.code.evaluate_functions(...
                obj.routines.planner_loss_commitment_discount,...
                ssdata(:,s0),xss,ss(:,s0),params(:,s0),def{s0},s0,s0);
            
            if ~utils.error.valid(lcd)
                
                retcode=6;
                
            end
            
            if ~retcode
                
                planner.objective{s0}=lcd(1);
                
                planner.commitment{s0}=lcd(2);
                
                planner.discount{s0}=lcd(3);
                
                if obj.is_optimal_simple_rule_model
                    
                    ww=zeros(obj.routines.planner_osr_support.size);
                    
                    ww(obj.routines.planner_osr_support.map)=lcd(4:end);
                    
                    ww=ww(obj.routines.planner_osr_support.partitions);
                    
                    ww=reshape(ww,orig_endo_nbr,orig_endo_nbr);
                    
                    good_order=obj.routines.planner_osr_support.derivatives_re_order;
                    
                    planner.weights{s0}=sparse(ww(good_order,good_order));
                    
                end
                
            end
            
        end
        
    end
    
    structural_matrices.planner=planner;
    
    if obj.is_optimal_simple_rule_model
        % change the ordering of the weight for the user. OSR will use
        % the structural matrices, which are ordered according to
        % order_var
        iov=obj.inv_order_var;
        
        for s0=1:h
            
            planner.weights{s0}=planner.weights{s0}(iov,iov);
            
        end
        
    end
    
    obj.solution.planner=planner; clear planner
    
end

    function [ys,nind]=forms_used_in_computation_of_derivatives()
        
        % spit out the forms to be used for the computation of derivatives
        %-----------------------------------------------------------------
        [the_leads,the_lags,nind]=...
            dsge_tools.create_endogenous_variables_indices(obj.lead_lag_incidence.before_solve);
        
        ys=zeros(nind,h);
        % derivatives taken wrt y+|y0|y-|shocks
        log_deriv_coefs=ones(nind+nx,h);
        
        is_log_var=obj.endogenous.is_log_var;
        
        % order of differentiation is different from alphabetic order
        %-------------------------------------------------------------
        yindex_deriv_order=obj.lead_lag_incidence.before_solve(obj.order_var,:);
        
        yindex_deriv_order=nonzeros(yindex_deriv_order(:))';
        
        long_is_log_var=is_log_var;
        
        long_is_log_var=[long_is_log_var(the_leads),long_is_log_var,long_is_log_var(the_lags)];
        
        for s11=1:h
            
            bgp=obj.solution.bgp{s11};
            
            % taking the approximation around the data and not the steady
            % state !
            sscurr=ssdata(:,s11);
            
            sslead=balanced_growth_path_powers(ssdata(:,s11),1);
            
            sslag=balanced_growth_path_powers(ssdata(:,s11),-1);
            
            ys(:,s11)=[sslead(the_leads);sscurr;sslag(the_lags)];
            
            tmp=ys(:,s11);
            
            tmp(~long_is_log_var)=1;
            
            % reorder according to the differentiation order
            log_deriv_coefs(1:nind,s11)=tmp(yindex_deriv_order);
            
        end
        
        structural_matrices.log_deriv_coefs=log_deriv_coefs;
        
        function [sstime]=balanced_growth_path_powers(sstime,c)
            
            sstime(is_log_var)=sstime(is_log_var).*bgp(is_log_var).^c;
            
            sstime(~is_log_var)=sstime(~is_log_var)+c*bgp(~is_log_var);
            
        end
        
    end

end