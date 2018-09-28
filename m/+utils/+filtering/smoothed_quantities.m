function [smoothed,shocks,pai0,ordered_endos,start_date]=smoothed_quantities(m,filtering)

ov=m.order_var;

ordered_endos=m.endogenous.name(ov);

sv=ts.collect(filtering.Expected_smoothed_variables);

start_date=serial2date(sv.start_date_number);

smoothed=double(sv(ordered_endos)).';

shks=ts.collect(filtering.Expected_smoothed_shocks);

shocks=double(shks).';

probs=ts.collect(filtering.smoothed_regime_probabilities);

pai0=double(probs);

pai0=pai0(1,:).';

end
