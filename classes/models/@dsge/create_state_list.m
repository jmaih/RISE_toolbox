function final_list=create_state_list(m,orders)
% create_state_list creates the list of the state variables in the solution
%
% Syntax
% -------
% ::
%
%   final_list=create_state_list(m)
%   final_list=create_state_list(m,orders)
%
% Inputs
% -------
%
% - **m** [dsge|rise] : model object
%
% - **orders** [integer array|{1:m.options.solve_order}] : approximation
%   orders
%
% Outputs
% --------
%
% - **final_list** [cellstr] : list of the state variables
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if isempty(m)
    final_list=struct();
    return
end

if nargin<2
    orders=1:m.options.solve_order;
end
orders=sort(orders);
shock_horizon=max(m.exogenous.shock_horizon);
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
final_list={};
old_state={};
for io=1:max(orders)
    new_state={};
    if isempty(old_state)
        new_state=state_list;
    else
        for istate=1:numel(old_state)
            new_state=[new_state,strcat(old_state{istate},',',state_list)];
        end
    end
    old_state=new_state;
    if any(io-orders==0)
        % save only the requested orders
        final_list=[final_list,new_state];
    end
end
end

