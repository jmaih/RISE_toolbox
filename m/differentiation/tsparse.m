classdef tsparse
    % tsparse -- constructs the transpose of a sparse matrix so as to save
    % memory
    properties
        v
        data
    end
    methods
        function obj=tsparse(ii,jj,vv,nrows,ncols)
            % tsparse -- constructor for tsparse objects
            %
            % Syntax
            % -------
            % ::
            %   obj=tsparse(ii,jj,vv,nrows,ncols)
            %
            % Inputs
            % -------
            % 
            % - **ii** [vector]: index for rows
            %
            % - **jj** [vector]: index for columns
            %
            % - **vv** [vector]: values
            %
            % - **nrows** [integer]: number of rows
            %
            % - **ncols** [integer]: number of columns
            %
            % Ouputs
            % -------
            %
            % - **obj** [tsparse]: tsparse object containing the transpose
            % of the matrix of interest
            %
            % More About
            % ------------
            %
            % Examples
            % ---------
            %
            % See also:
            if nargin
                obj.data=sparse(jj,ii,vv,ncols,nrows);
                obj.v=vv;
            end
        end
        function c=mtimes(a,b)
            a_tsp=isa(a,'tsparse');
            b_tsp=isa(b,'tsparse');
            if a_tsp && ~b_tsp
                c=a.data.'*b;
            elseif ~a_tsp && b_tsp
                c=a*b.data.';
            else
                c=a.data.'*b.data.';
            end
        end
        function varargout=size(A,d)
            nout=nargout();
            varargout=cell(1,nout);
            if nargin==1
                for iout=1:nout
                    if iout==1
                        varargout{iout}=size(A.data,2);
                    elseif iout==2
                        varargout{iout}=size(A.data,1);
                    else
                        varargout{iout}=1;
                    end
                end
            else
                if nout~=1
                    error('with two input arguments, the number of output arguments can only be 1')
                end
                w=(d==1)*2+(d==2)*1;
                if w==0
                    varargout{1}=1;
                else
                    varargout{1}=size(A.data,w);
                end
            end
        end
        function d=double(a)
            d=a.data;
        end
        function C=reshape(A,r,c)
            % the function produces reshape(transpose(A.data),r,c)
            tvec=@(x)vec(transpose(x));
            A=A.data;
            [ra,ca]=size(A);
            Af=[];
            ss=memory;
            max_bytes=ss.MaxPossibleArrayBytes/4;
            % MemAvailableAllArrays = MaxPossibleArrayBytes so we divide by
            % 4 above to be economical...
            ncell=200;
            C=cell(200,1);
            iter=0;
            ic=0;
            while ic<c
                % choose the number of rows
                %--------------------------
                nrows=ceil(max_bytes/ca);
                nrows=min(ra,nrows);
                
                swallow_rows();
                nf=numel(Af);
                if nf<r
                    error('reshaping fails: probably too large matrix')
                end
                nc=floor(nf/r);
                r_nc=r*nc;
                [ii,jj,vv]=find(reshape(Af(1:r_nc),r,nc));
                if ~isempty(ii)
                    iter=iter+1;
                    if iter==ncell
                        % increase number of cells
                        ncell=ncell+200;
                        C(ncell)={};
                    end
                    C{iter}=[ii(:),jj(:)+ic,vv(:)];
                end
%                 C(:,ic+(1:nc))=reshape(Af(1:r_nc),r,nc);
                ic=ic+nc;
                Af=Af(r_nc+1:end);
                A(1:nrows,:)=[];
                ra=ra-nrows;
            end
            C=cell2mat(C(1:iter));
            C=sparse(C(:,1),C(:,2),C(:,3),r,c);
            function swallow_rows()
                Af=[Af;tvec(A(1:nrows,:))];
            end
        end
        function B=subsref(A,S)
            switch S.type
                case '()'
                    ncols=numel(S.subs);
                    S.subs=S.subs(1:min(2,ncols));
                    S.subs=fliplr(S.subs);
                    B=subsref(A.data,S);
                    B=transpose(B);
                case '.'
                    B=A.(S.subs);
                case '{}'
                    error('only subsref of types () and . are supported')
            end
        end
    end
end