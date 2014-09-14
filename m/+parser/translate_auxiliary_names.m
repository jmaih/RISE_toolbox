function out=translate_auxiliary_names(names)
out=regexprep(names,'(\w+)_AUX_L_(\d+)','$1{-$2}');
out=regexprep(out,'(\w+)_AUX_F_(\d+)','$1{+$2}');
end