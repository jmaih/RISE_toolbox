% PAD Insert leading and trailing spaces
%    NEWSTR = PAD(STR) inserts space characters at the end of each element
%    of STR to match the width of the longest string. STR can be a string
%    array, character vector, or cell array of character vectors.
% 
%    NEWSTR = PAD(STR,WIDTH) inserts space characters to a width specified
%    by WIDTH. Elements in STR that are longer than WIDTH are unchanged. The
%    default WIDTH is max(strlength(STR)).
% 
%    NEWSTR = PAD(STR,SIDE) inserts space characters to STR on the specified
%    SIDE to match the width of the longest string. SIDE can be 'left',
%    'right', or 'both'. The default SIDE is 'right'.
% 
%    NEWSTR = PAD(STR,PAD_CHARACTER) inserts PAD_CHARACTER at the end of
%    each string element in STR to match the width of the longest string.
%    PAD_CHARACTER must be exactly one character. The default PAD_CHARACTER
%    is a space character, ' '.
% 
%    NEWSTR = PAD(STR,WIDTH,SIDE) inserts space characters to the specified
%    SIDE and WIDTH.
% 
%    NEWSTR = PAD(STR,WIDTH,PAD_CHARACTER) inserts PAD_CHARACTER to the
%    specified WIDTH.
% 
%    NEWSTR = PAD(STR,WIDTH,SIDE,PAD_CHARACTER) inserts PAD_CHARACTER to the
%    specified SIDE and WIDTH.
% 
%    Example:
%        STR = string({'head';'shoulders';'knees';'toes'});
%        pad(STR)
% 
%        returns
% 
%            "head     "
%            "shoulders"
%            "knees    "
%            "toes     "
%    
%    Example:
%        STR = string({'head';'shoulders';'knees';'toes'});
%        pad(STR,10,'left','.')
% 
%        returns
% 
%            "......head"
%            ".shoulders"
%            ".....knees"
%            "......toes"
%        
%    See also STRIP, STRING/PLUS
%
%    Reference page in Doc Center
%       doc pad
%
%    Other functions named pad
%
%       codistributed/pad    tall/pad
%