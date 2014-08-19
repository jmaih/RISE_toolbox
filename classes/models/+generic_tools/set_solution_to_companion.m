function [A,B]=set_solution_to_companion(obj)
endo_nbr=obj.endogenous.number(end);
exo_nbr=sum(obj.exogenous.number);
reg_nbr=obj.markov_chains.regimes_number;
if isa(obj,'svar')
    A=obj.solution.m_x;
    B=obj.solution.m_e;
    cc=(obj.nlags-1)*endo_nbr;
    for ireg=1:reg_nbr
        A{ireg}=[A{ireg}(:,1:obj.nlags*endo_nbr);
            eye(cc),zeros(cc,endo_nbr)];
        B{ireg}=[B{ireg};
            zeros(cc,exo_nbr)];
    end
elseif isa(obj,'dsge')
    [A,B]=inv_order_var_solution();
end
    function [Tz,Re]=inv_order_var_solution()
        %         ov=obj.order_var.after_solve;
        iov=obj.inv_order_var.after_solve;
        z_pb=obj.locations.after_solve.z.pb;
        t_pb=obj.locations.after_solve.t.pb;
        e_0=obj.locations.after_solve.z.e_0;
        Re=cell(1,reg_nbr);
        Tz=Re;
        tmp=zeros(endo_nbr);
        for isol=1:reg_nbr
            tmp(:,t_pb)=obj.solution.Tz{isol}(:,z_pb);
            % separate autoregressive part from shocks
            % the rows are already ordered in the iov order and so we
            % re-order the columns only
            %--------------------------------------------------------------
            Tz{isol}=tmp(:,iov);
            Re{isol}=obj.solution.Tz{isol}(:,e_0(1):end);
            if isol==1
                npges=size(Re{isol},2)/exo_nbr;
            end
            Re{isol}=reshape(Re{isol},[endo_nbr,exo_nbr,npges]);
        end
    end

end
