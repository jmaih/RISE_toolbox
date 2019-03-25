% CHAR Create character array
%    CHR = CHAR(X) converts array X of nonnegative integer codes into a
%    character array. Valid codes range from 0 to 65535, where codes 0 to
%    127 correspond to 7-bit ASCII characters. The characters that MATLAB
%    can process (other than 7-bit ASCII characters) depend on your current
%    locale setting. Use DOUBLE to convert characters to numeric codes.
% 
%    CHR = CHAR(C), when C is a cell array of character vectors, places each 
%    element of C into a row of the character array CHR. Use CELLSTR to
%    convert back.
% 
%    CHR = CHAR(STR), when STR is a string array, converts each element of STR into 
%    a row of the character array CHR. Use STRING to convert back.
% 
%    CHR = CHAR(T1,T2,T3,..) forms the character array CHR containing the
%    text from T1,T2,T3,... as rows. CHAR automatically pads each row with
%    spaces in order to form a character array. Each text parameter, Ti,
%    can itself be a character array. This allows the creation of
%    arbitrarily large character arrays. If Ti has no characters, then the
%    corresponding row of CHR is filled with spaces.
% 
%    See also STRING, DOUBLE, CELLSTR, ISCELLSTR, ISCHAR, ISSTRING.
%
%    Reference page in Doc Center
%       doc char
%
%    Other functions named char
%
%       categorical/char      inline/char    splanar/char    tall/char
%       codistributed/char    opaque/char    sym/char
%