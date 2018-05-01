function output=append_file(output,rline,definitions)
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


rf=include_file();
output=[output;rf];
    function rf=include_file()
        quotes=strfind(rline.rest,'"');
        if isempty(rline.rest)
            error([mfilename,':: missing quotes after statement @#include in ',...
                rline.FileName,' at line ',...
                rline.row_number_string])
        end
        if numel(quotes)~=2
            error([mfilename,':: expecting a pair of quotes in ',...
                rline.FileName,' at line ',...
                rline.row_number_string])
        end
        newfile=strtrim(rline.rest(quotes(1)+1:quotes(2)-1));
        if ~exist(newfile,'file')
            error([mfilename,':: file ',newfile,' not found::  ',...
                rline.FileName,...
                ' at line ',rline.row_number_string])
        end
        rf=parser.preparse(newfile,definitions);
    end
end
