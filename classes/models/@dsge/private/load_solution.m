function [T,R,steady_state,state_vars_location]=load_solution(obj,type)
% when the type is ov, R is empty while T is of size solve_order x
% regimes_number. Basically, T is returned as Tz, Tzz, Tzzz, etc. where the
% rows that were originally in alphabetical order are put in order of
% simulation/processing
% when the type is ov, R corresponds to the shocks and T corresponds to the
% autoregressive elements and the matrices in T are square

if ~any(strcmp(type,{'ov','iov'}))
    error('type should be ov (order var) or iov (inv order var)')
end
is_alphabetical_order=strcmp(type,'iov');

regimes_number=obj.markov_chains.regimes_number;
ov=obj.order_var.after_solve;
order=obj.options.solve_order;
T=cell(order,regimes_number);
R=cell(1,regimes_number);
steady_state=obj.solution.ss;
state_vars_location=obj.locations.after_solve.t.pb;
% update T and steady state
%--------------------------
order_var_solution();

% now get the inv_order_var solution if necessary
%------------------------------------------------
if is_alphabetical_order
    iov=obj.inv_order_var.after_solve;
    % only for order 1
    inv_order_var_solution()
end
    function inv_order_var_solution()
        endo_nbr_=obj.endogenous.number(end);
        exo_nbr_=sum(obj.exogenous.number);
        z_pb=obj.locations.after_solve.z.pb;
        t_pb=state_vars_location;
        state_vars_location=[];
        e_0=obj.locations.after_solve.z.e_0;
        Tz=R;
        tmp=zeros(endo_nbr_);
        for isol=1:regimes_number
            % this solution comes from order_var below and so we re-order
            % both rows and columns
            tmp(:,t_pb)=T{1,isol}(:,z_pb);
            % separate autoregressive part from shocks
            %-----------------------------------------
            Tz{isol}=tmp(iov,iov); 
            R{isol}=T{1,isol}(iov,e_0(1):end);
            if regimes_number>1
                if isol==1
                    npges=size(R{isol},2)/exo_nbr_;
                end
                R{isol}=reshape(R{isol},[endo_nbr_,exo_nbr_,npges]);
            end
        end
        T=Tz;
        if order>1
            warning([mfilename,':: Only first order will be used'])
        end
    end

    function order_var_solution()
        zzz=repmat('z',1,order);
        for io=1:order
            for ireg=1:regimes_number
                if ~is_alphabetical_order && io==1
                    % do this only if the alphabetical order is not needed
                    steady_state{ireg}=steady_state{ireg}(ov);
                end
                % re-order the rows
                %------------------
                T{io,ireg}=obj.solution.(['T',zzz(1:io)]){ireg}(ov,:);
            end
        end
    end
end