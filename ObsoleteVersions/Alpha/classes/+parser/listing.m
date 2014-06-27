function xout=listing()

xout=struct('name',{},'tex_name',{},... all
    'is_in_use',{},... exogenous and parameters
    'governing_chain',{},'is_switching',{},'is_measurement_error',{},'is_trans_prob',{},... parameters
    'is_log_var',{},... endogenous
    'is_endogenous',{},... observables
    'max_lead',{},... variables and parameters
    'max_lag',{},... variables and parameters
    'state_id',{},... observables
	'is_auxiliary',{}... endogenous
    );
