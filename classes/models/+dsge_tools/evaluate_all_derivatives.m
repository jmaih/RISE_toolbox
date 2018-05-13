function [structural_matrices,retcode]=evaluate_all_derivatives(obj,structural_matrices,ssdata)
% evaluate all the derivatives up to the desired order: analytical
% derivatives can be evaluated sequentially, but for algorithmic
% derivatives and for numerical derivatives, it is more economical to do
% them at once, instead of calling the same functions several times

if nargin < 3
    
    ssdata = [];
    
end

solve_time_consistent_switch=obj.options.solve_time_consistent_switch;

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
    
    if obj.is_optimal_policy_model
        % loose commitment may call for this
        %-----------------------------------        
        ssdata=ssdata(:,1);
        
    end
    
    if size(ssdata,2)~=1
        
        error('wrong number of columns for ssdata')
        
    else
        
        ssdata=ssdata(:,ones(h,1));
        
    end
    
end

if size(ssdata,1)~=size(ss,1)
    
    error('wrong number of rows for ssdata')
    
end

[ys,nind,lead_pos]=forms_used_in_computation_of_derivatives();

if is_has_data
    % update residuals with the data
    for ireg=1:h
        
        s0=ireg;
        
        s1=ireg;
        
        ys01=yvector(s0,s1);
        
        structural_matrices.user_resids(:,ireg)=utils.code.evaluate_functions(...
            obj.routines.probs_times_dynamic,...
            ys01,xss,ss(:,s0),params(:,s0),def{s0},s0,s1);
        
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
        
         ys01=yvector(s0,s1);
        
        if ~retcode
            % Note: G(s0,s1) =: ps0(s0,s1)*F(s0)
            if symbolic_type
                
                if s1==1 && s0==1
                    
                    tmp_=obj.routines.probs_times_dynamic_derivatives;
                    
                    if isa(tmp_,'function_handle')
                        
                        max_order=nargout(tmp_);
                        
                    else
                        
                        max_order=numel(tmp_);
                        
                    end
                    
                    clear tmp_
                    
                    if solve_order>max_order
                        
                        error(['Perturbation of order ',int2str(solve_order),...
                            ' requested but symbolic derivatives available ',...
                            'only up to order ',int2str(max_order),...
                            '. Compute higher-order symbolic derivatives or ',...
                            'switch to automatic derivatives'])
                        
                    end
                    
                end
                
                [G01{1:solve_order}]=utils.code.evaluate_functions(obj.routines.probs_times_dynamic_derivatives,...
                    ys01,xss,ss(:,s0),params(:,s0),def{s0},s0,s1);
                
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
                    ys01,xss,ss(:,s0),params(:,s0),def{s0},s0,s1);
                
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
                    ys01,xss,ss(:,s0),params(:,s0),def{s0},s0,s1);
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
                        multiply_bsxfun(G01{io},zkz);
                    
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
    
    endo_nbr=obj.endogenous.number;
    
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
                
                planner.weights{s0}=calculate_weights(lcd,s0);
                
            end
            
        end
        
    end
    
    structural_matrices.planner=planner;
    
end

    function y01=yvector(s0,s1)
        
        y01=ys(:,s0);
        
        if solve_time_consistent_switch
            
            y01(lead_pos)=ys(lead_pos,s1);
            
        end
        
    end

    function w=calculate_weights(lcd,s0)
        
        if symbolic_type
            
            ww=zeros(obj.routines.planner_osr_support.size);
            
            ww(obj.routines.planner_osr_support.map)=lcd(4:end);
            
            ww=ww(obj.routines.planner_osr_support.partitions);
            
            nwrt=obj.routines.planner_osr_support.nwrt;
            
            ww=reshape(ww,nwrt,nwrt);
            
            w=zeros(endo_nbr);
            
            good_order=obj.routines.planner_osr_support.derivatives_re_order;
            
            w(good_order,good_order)=ww;
        
        elseif automatic_type
            
            w=automatic_policy_weigths(obj,...
                ssdata(:,s0),xss,ss(:,s0),params(:,s0),def{s0},s0,s0);
            
        else
            
            iov=obj.inv_order_var;
            
            w(iov,iov)=utils.numdiff.hessian(obj.routines.planner_objective{1},...
                ssdata(:,s0),[],xss,ss(:,s0),params(:,s0),def{s0},s0,s0);
            
        end
        
        w=sparse(w);
        
    end

    function [ys,nind,lead_positions,lag_positions]=...
            forms_used_in_computation_of_derivatives()
        
        % spit out the forms to be used for the computation of derivatives
        %-----------------------------------------------------------------
        [the_leads,the_lags,nind]=...
            dsge_tools.create_endogenous_variables_indices(obj.lead_lag_incidence.before_solve);
        
        lead_positions=1:numel(the_leads);
        
        curr=1:size(ssdata,1);
        
        lag_positions=numel(the_leads)+numel(curr)+(1:numel(the_lags));
        
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

function G=multiply_bsxfun(G,zkz)

if all(zkz==1)
    
    return
    
end

try
    
    G=bsxfun(@times,G,zkz);
    
catch
    
    for ii=1:size(G,1)
        
        tmp=find(G(ii,:));
        
        G(ii,tmp)=G(ii,tmp).*zkz(tmp);
        
    end
    
end

end

function W=automatic_policy_weigths(obj,y,xss,ss,params,def,s0,s1) %#ok<INUSL,INUSD>

sm=func2str(obj.routines.planner_objective{1});

rp=find(sm==')',1,'first');

sm=sm(rp+1:end);

sm=parser.symbolic_model(sm);

tokens=regexp(sm,'\<(((y|param|ss)_\d+)|s0|s1)\>','match');

f=parser.cell2matize(tokens);

f=strrep(f,'|',',');

f=str2func(['@',f,sm]);

nt=numel(tokens);

pos=nan(nt,1);

types=nan(nt,1);

isactive=false(1,nt);

for ii=1:nt
    
    tok=tokens{ii};
    
    if strcmp(tok(1),'y')
        
        pos(ii)=str2double(tok(3:end));
        
        types(ii)=1;
        
        isactive(ii)=true;
        
    elseif strcmp(tok(1),'p')
        
        pos(ii)=str2double(tok(7:end));
        
        types(ii)=2;
        
    elseif strcmp(tok(1),'s')
        
        pos(ii)=str2double(tok(4:end));
        
        types(ii)=3;
        
    elseif strcmp(tok(1),'d')
        
        pos(ii)=str2double(tok(5:end));
        
        types(ii)=4;
        
    else
        
        disp([tok,': unrecognized'])
        
    end
    
end

ping=pos(types==1);

n=numel(ping);

% push values
%-------------
vals=nan(nt,1);

vals(types==1)=y(ping);

vals(types==2)=params(pos(types==2),1);

vals(types==3)=ss(pos(types==3));

vals(types==4)=def(pos(types==4));

toks=[tokens(:),num2cell(vals(:))];

% Take automatic derivatives
%----------------------------
der=aplanar.diff({f},toks(isactive,:),toks(~isactive,:),2);

W=zeros(obj.endogenous.number);

W(ping,ping)=reshape(full(der{2}),n,n);

ov=obj.order_var;

W=W(ov,ov);

end