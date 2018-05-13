function [final_list,kept]=create_state_list(m,orders,compact_form)
% create_state_list creates the list of the state variables in the solution
%
% ::
%
%
%   final_list=create_state_list(m)
%   final_list=create_state_list(m,orders)
%
% Args:
%
%    - **m** [dsge|rise] : model object
%
%    - **orders** [integer array|{1:m.options.solve_order}] : approximation
%      orders
%
%    - **compact_form** [true|{false}] : if true, only unique combinations
%      will be returned. Else, all combinations will be returned.
%
% Returns:
%    :
%
%    - **final_list** [cellstr] : list of the state variables
%
%    - **kept** [vector] : location of kept state variables (computed only if
%      compact_form is set to true)
%
% Note:
%
% Example:
%
%    See also:

if isempty(m)
    
    final_list=cell(0,4);
    
    return
    
end

if nargin<3
    
    compact_form=false;
    
    if nargin<2
    
        orders=[];
    
    end
    
end

if isempty(orders)

    orders=1:m.options.solve_order;

end

orders=sort(orders);

shock_horizon=max(m.exogenous.shock_horizon(:));

exo_list=get(m,'exo_list');

% "predetermined" and then "both" variables
%-------------------------------------------
state_list = [get(m,'endo_list(predetermined)'),get(m,'endo_list(pred_frwrd_looking)')];

% lag those names
%-----------------
state_list = parser.lag_names(state_list); 

% add the perturbation parameter then the exogenous
%--------------------------------------------------
state_list=[state_list,'@sig',exo_list];

% finally the future shocks
%--------------------------
for ik=1:shock_horizon

    state_list=[state_list,strcat(exo_list,'{+',int2str(ik),'}')];

end

n=numel(state_list);

final_list={};

old_state={};

kept=cell(numel(orders),1);

if compact_form
    
    kept=utils.kronecker.shrink_expand(n,max(orders));
    
end

for io=1:max(orders)

    new_state={};
    
    if isempty(old_state)
    
        new_state=state_list;
    
    else
        
        for istate=1:numel(old_state)
        
            new_state=[new_state,strcat(old_state{istate},',',state_list)]; %#ok<*AGROW>
        
        end
        
    end
    
    old_state=new_state;
        % save only the requested orders
    if any(io-orders==0)
   
        if compact_form
        
            new_state=new_state(kept{io});
        
        end
        
        final_list=[final_list,new_state];
    
    end
    
end

kept=cell2mat(kept(orders));

end