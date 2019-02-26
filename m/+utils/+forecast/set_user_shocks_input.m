function user_input=set_user_shocks_input(exo_names,simul_user_shocks_draws)

% e.g. fn=@(n)(sum(randn(k,n).^2,1)-k)/sqrt(2*k);

% e.g. xlist.E_cons_pref='cons_pref_draw'
% e.g. xlist.E_cons_pref={'cons_pref_draw',1,2,3,4,5}
% e.g. xlist.E_wage_markup={@(n,k)(sum(randn(k,n).^2,1)-k)/sqrt(2*k),5}
% e.g. xlist.E_price_markup=@(n)(sum(randn(3,n).^2,1)-3)/sqrt(2*3)

if isempty(simul_user_shocks_draws)
    
    user_input=[];
    
    return
    
end

user_list=fieldnames(simul_user_shocks_draws);

locs=locate_variables(user_list,exo_names);

n=numel(user_list);

user_input=cell(n,2);

for ii=1:n
    
    shk=user_list{ii};
    
    user_input(ii,:)={locs(ii),simul_user_shocks_draws.(shk){2}};
    
end

end