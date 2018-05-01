function xout=string_mult(a,b)
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


if ~ischar(a)
    a=char(a);
end
if ~ischar(b)
    b=char(b);
end
a(isspace(a))=[];
b(isspace(b))=[];
% check whether there is a semicolon at the end of any of the strings
%--------------------------------------------------------------------
is_semicol_a=strcmp(a(end),';');
is_semicol_b=strcmp(b(end),';');
is_semicol=is_semicol_a||is_semicol_b;
if is_semicol_a
    a(end)=[];
end
if is_semicol_b
    b(end)=[];
end

if strcmp(a,'0')||strcmp(b,'0')
    xout='0';
elseif strcmp(a,'1')
    xout=b;
elseif strcmp(a,'1')
    xout=b;
else
    xout=['(',a,')*(',b,')'];
end
if is_semicol
    xout=[xout,';'];
end
end