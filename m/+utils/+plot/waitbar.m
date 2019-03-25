% WAITBAR Display wait bar.
%    H = WAITBAR(X,'message', property, value, property, value, ...)
%    creates and displays a waitbar of fractional length X.  The
%    handle to the waitbar figure is returned in H.
%    X should be between 0 and 1.  Optional arguments property and
%    value allow to set corresponding waitbar figure properties.
%    Property can also be an action keyword 'CreateCancelBtn', in
%    which case a cancel button will be added to the figure, and
%    the passed value string will be executed upon clicking on the
%    cancel button or the close figure button.
% 
%    WAITBAR(X) will set the length of the current bar to the fractional
%    length X.
% 
%    WAITBAR(X,H) will set the length of the bar in waitbar H
%    to the fractional length X.
% 
%    WAITBAR(X,H,'message') will update the message text in
%    the waitbar figure, in addition to setting the fractional
%    length to X.
% 
%    WAITBAR is typically used inside a FOR loop that performs a
%    lengthy computation.
% 
%    Example:
%        h = waitbar(0,'Please wait...');
%        for i=1:1000,
%            % computation here %
%            waitbar(i/1000,h)
%        end
% 
%    See also DIALOG, MSGBOX.
%
%    Reference page in Doc Center
%       doc waitbar
%
%