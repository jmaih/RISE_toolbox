function rawline_=remove_comments(rawline_)
% remove_comments - removes comments of type "//", "%"
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% - "\%" is not considered a comment. This is the way to mask "%"
%
% Examples
% ---------
%
% See also: 

% locate comments of type %: those are the ones not preceeded with a \
%----------------------------------------------------------------------
loc_=regexp(rawline_,'(?<!\\)%'); % strfind(rawline_,'%');

if ~isempty(loc_)
    
    rawline_=rawline_(1:loc_(1)-1);
    
end

loc_=strfind(rawline_,'//');

if ~isempty(loc_)
    
    rawline_=rawline_(1:loc_(1)-1);
    
end

end
