%--- help for contains ---
%
% CONTAINS True if pattern is found in text.
%    TF = CONTAINS(STR,PATTERN) returns 1 (true) if STR contains PATTERN,
%    and returns 0 (false) otherwise.
% 
%    STR can be a string array, a character vector, or a cell array of
%    character vectors. So can PATTERN. PATTERN and STR need not be the same
%    size. If PATTERN is a string array or cell array, then CONTAINS returns
%    true if it finds any element of PATTERN in STR. If STR is a string
%    array or cell array, then TF is a logical array that is the same size.
%  
%    TF = CONTAINS(STR,PATTERN,'IgnoreCase',IGNORE) ignores case when searching 
%    for PATTERN in STR if IGNORE is true. The default value of IGNORE is false.
%  
%    Examples
%        STR = "data.tar.gz";
%        P = "tar";
%        contains(STR,P)                   returns  1
% 
%        STR = ["abstracts.docx","data.tar.gz"];
%        P = 'tar';         
%        contains(STR,P)                   returns  [0 1]
% 
%        STR = 'data.tar.gz';
%        P = {'docx','tar'};
%        contains(STR,P)                   returns  1
% 
%        STR ={'DATA.TAR.GZ','SUMMARY.PPT'};
%        P = "tar";
%        contains(STR,P,'IgnoreCase',true) returns  [1 0]
% 
%    See also endsWith, startsWith.
%
%    Reference page in Doc Center
%       doc contains
%
%