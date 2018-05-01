function xout=kron(a,b)
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

% kronecker multiplication of cell arrays of strings

if ischar(a)
    a=cellstr(a);
end
if ischar(b)
    b=cellstr(b);
end
[ra,ca]=size(a);
[rb,cb]=size(b);
xout=cell(ra*rb,ca*cb);
proto_b=cell(rb,cb);
for irow_a=1:ra
    stretch_ra=(irow_a-1)*rb+(1:rb);
    for icol_a=1:ca
        arc=a{irow_a,icol_a};
        stretch_ca=(icol_a-1)*cb+(1:cb);
        for irow_b=1:rb
            for icol_b=1:cb
                proto_b{irow_b,icol_b}=parser.string_mult(arc,b{irow_b,icol_b});
            end
        end
        xout(stretch_ra,stretch_ca)=proto_b;
    end
end
end
