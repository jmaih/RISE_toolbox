function out=translate_auxiliary_names(names)

% variables with lags
out=regexprep(names,'(\w+)_AUX_L_(\d+)','$1{-$2}');

% variables with leads
out=regexprep(out,'(\w+)_AUX_F_(\d+)','$1{+$2}');

end