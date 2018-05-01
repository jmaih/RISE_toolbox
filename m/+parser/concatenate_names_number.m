function out=concatenate_names_number(names,lag)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

test=@myadd; %#ok<NASGU>
% translate auxiliaries
%----------------------
 out=parser.translate_auxiliary_names(names);
 
% set the lead/lag to 0 for variables that are contemporaneous
%-------------------------------------------------------------
out=regexprep(out,'(\w+)(?!({|}|\w+))','$1{0}');

% update the lead/lag
%--------------------
out=regexprep(out,'(\w+){((\+|-)?\d+)}','$1{${test($2)}}');

% remove any lead/lag equal to 0
%-------------------------------
out=strrep(out,'{0}','');

% add a plus sign for leads
%--------------------------
out=regexprep(out,'(\w+){(\d+)}','$1{+$2}');

    function out=myadd(x)
        out=sprintf('%0.0f',str2double(x)+lag);
    end

end