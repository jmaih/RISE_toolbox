function status=pdflatex(filename)
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

% get the compiler from the rise_root
compiler=getappdata(0,'rise_pdflatex');
if exist([filename,'.pdf'],'file')
    delete([filename,'.pdf'])
end
the_string=[compiler,' ',filename];
% close all hanging fids
fclose('all');
retcode=system(the_string);
if ~retcode
    % run again in order to make sure the references are shown
    system(the_string);
    system(the_string);
end
if nargout
    status=retcode;
end

end