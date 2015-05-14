function [mompos,sigpos]=moments_positions(siz,order)
% moments_positions -- finds the coordinates of moments in Ekron(u,u,...,u)
%
% Syntax
% -------
% ::
%
%   [mompos,sigpos]=moments_positions(siz,order)
%
% Inputs
% -------
%
% - siz : [struct] with fields
%   - np: number of predetermined
%   - nb: number of predetermined and forward-looking
%   - ne: number of shocks
%   - nz: (total) number of state variables (pred,bobth,sig,shocks)
%
% - order : [numeric] order of the moments
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

shocks_id=siz.np+siz.nb+1+(1:siz.ne);
mompos=engine(shocks_id);
sig_id=siz.np+siz.nb+1;
sigpos=engine(sig_id);

    function A=engine(locs)
        parts=siz.nz*ones(1,order);
        if order==1
            parts=[parts,1];
        end
        A=zeros(parts);
        for ii=locs
            if order==1
                A(ii)=1;
            else
                for jj=locs
                    if jj==ii
                        if order==2
                            A(ii,jj)=1;
                        else
                            for kk=locs
                                if kk==jj
                                    if order==3
                                        A(ii,jj,kk)=1;
                                    else
                                        for ll=locs
                                            if ll==kk
                                                if order==4
                                                    A(ii,jj,kk,ll)=1;
                                                else
                                                    for mm=locs
                                                        if mm==ll
                                                            if order==5
                                                                A(ii,jj,kk,ll,mm)=1;
                                                            else
                                                                error('approximation order greater than 5 not supported')
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        A=logical(vec(A));
    end
end
