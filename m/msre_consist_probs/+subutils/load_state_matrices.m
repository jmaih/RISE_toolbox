function [mat,Resids,Q,c,own_stuff]=load_state_matrices(m,structural_matrices,T,x0,inner_probabilities,control_shocks,own_stuff)

if nargin<6
    % this is to avoid using persistent variables
    own_stuff=[];
    
end

if isempty(own_stuff)
    
    own_stuff=struct();
    
    own_stuff.is_constant = isempty(inner_probabilities);
    
    own_stuff.derivs=struct();
    
end

if own_stuff.is_constant
    
    c=[];
    
    if ~isempty(own_stuff.derivs.A00)
        
        mat.Aplus=own_stuff.derivs.Aplus0;
        
        mat.A0=own_stuff.derivs.A00;
        
        mat.Aminus=own_stuff.derivs.Aminus0;
        
        mat.B=own_stuff.derivs.B0;
        
        Resids=own_stuff.Resids;
        
        Q=own_stuff.Q;
        
        return
        
    end
    
else
    
    c=simulated_probabilities(T,x0,inner_probabilities,control_shocks);
    
    m=set(m,'parameters',[inner_probabilities(:,1),num2cell(c)']);
    
end

[structural_matrices,retcode]=dsge_tools.evaluate_all_derivatives(m,structural_matrices,[]);

[pos,siz,shock_horizon]=dsge_tools.rehash_topology(m,structural_matrices);

[sm]=utils.solve.pull_first_order_partitions(structural_matrices.dv,pos.v);

Resids=structural_matrices.user_resids;

adjusted=struct();

adjusted.bf_cols=pos.t.bf;

adjusted.pb_cols=pos.t.pb;

adjusted.nd=siz.nd;

adjusted.siz=siz; % adjusted sizes

accelerate=false;
% accelerate=options.solve_accelerate && siz.ns;

[sm,adjusted,accelSupport]=dsge_tools.aggregate_matrices(sm,siz,adjusted,accelerate); %#ok<ASGLU>

[mat.Aplus,mat.A0,mat.Aminus]=dsge_tools.full_state_matrices(siz,sm);

Q=update_transition_matrix_since_steady_state_is_not_resolved();

% Q = structural_matrices.transition_matrices.Q;

mat.B=full(sm.de_0{1});

mat.B=mat.B(:,:,ones(siz.h,1));

for istate = 2:siz.h
    
    mat.B(:,:,istate)=sm.de_0{istate};
    
end

if own_stuff.is_constant
    
    own_stuff.derivs.Aplus0=mat.Aplus;
    
    own_stuff.derivs.A00=mat.A0;
    
    own_stuff.derivs.Aminus0=mat.Aminus;
    
    own_stuff.derivs.B0=mat.B;
    
    own_stuff.Resids=Resids;
    
    own_stuff.Q=Q;
    
end

    function Q=update_transition_matrix_since_steady_state_is_not_resolved()
        
        y1=zeros(m.endogenous.number,1);
        
        p=m.parameter_values;
        
        d={nan};
        
        [TransMat,retcode]=compute_steady_state_transition_matrix(...
            m.routines.transition_matrix,y1,p(:,1),d{1},...
            sum(m.exogenous.number));
        
        Q=TransMat.Q;
        
        if retcode
            
            decipher(retcode)
            
        end
        
    end

end