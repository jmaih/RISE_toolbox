function add_str=lead_lag_string(lead_or_lag)
if lead_or_lag>0
    add_str='_AUX_F_';
elseif lead_or_lag<0
    add_str='_AUX_L_';
else
    error('input must be strictly positive or strictly negative')
end
end