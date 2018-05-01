function Z=A_times_B_kron_C(A,B,C)
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

% computes the A*kron(B,C)

[rb,cb]=size(B);
[rc,cc]=size(C);
[ra,ca]=size(A);
if ~ca==rb*rc
    error('wrong matrices sizes')
end
% compute A*kron(B,I)
 ABI=A_times_B_kron_I(true);

% compute D*kron(I,C)
Z=D_kron_I_C();

    function Z=D_kron_I_C()
        Z=zeros(ra,cb*cc);
        offset=0;
        for icol=1:cb
            Z(:,(icol-1)*cc+1:icol*cc)=ABI(:,offset+(1:rc))*C;
            offset=offset+rc;
        end
    end

    function ABI=A_times_B_kron_I(fast)
        % compute A*kron(B,I)
        if nargin==0
            fast=true;
        end
        if fast
            ABI=zeros(ra*rc,cb);
            for jcol=1:rc
                ABI((jcol-1)*ra+1:jcol*ra,:)=A(:,jcol:rc:end)*B;
            end
            ABI=reshape(ABI,ra,cb*rc);
        else
            ABI=zeros(ra,rb*rc);
            jj=0;
            for icol=1:cb
                offset=0;
                for jcol=1:rc
                    jj=jj+1;
                    ABI(:,jj)=A(:,offset+1:rc:end)*B(:,icol);
                    offset=offset+1;
                end
            end
        end
    end

end