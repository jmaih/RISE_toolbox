classdef rise_anonymous
    properties
        size
        indices
        func
    end
    methods
        function obj=rise_anonymous(siz,string,ind)
            if nargin
                obj.size=siz;
                obj.indices=ind;
                obj.func=string;
            end
        end
        function val=eval(obj,varargin)
            n=numel(obj);
            val=cell(1,n);
            for ii=1:n
                val{ii}=zeros(obj(ii).size);
                val{ii}(obj(ii).indices)=obj(ii).func(varargin{:});
            end
        end
        function [Q,retcode]=kron(obj,varargin)
            n=numel(obj);
            Q=eval(obj(1),varargin{:});
            Q=Q{1};
            for ii=2:n
                val_i=eval(obj(ii),varargin{:}); %#ok<*EVLC>
                val_i=val_i{1};
                if isequal(Q,1)
                    Q=val_i;
                elseif isequal(val_i,1)
                    % no need to update Q
                else
                    Q=kron(Q,val_i{1});
                end
            end
            retcode=0;
            if any(any(isnan(Q))) || any(any(Q<0)) || any(any(Q>1)) 
                Q=[];
                retcode=3;
% %             else
% %                 % protect against under/overflow
% %                 n=size(Q,1);
% %                 diag_terms=(0:n-1)*n+(1:n);
% %                 Q(diag_terms)=0;
% %                 Q(diag_terms)=1-sum(Q,2);
            end
            
        end
        function write2file(obj,filename)
            n=numel(obj);
            thedot=strfind(filename,'.');
            if ~isempty(thedot)
                filename=filename(1:thedot-1);
            end
            args=func2str(obj(1).func);
            right_par=find(args==')',1,'first');
            args=args(3:right_par-1);
            fid=fopen([filename,'.m'],'w');
            fprintf(fid,'%s\n',['function Q=',filename,'(',args,')']);
            fprintf(fid,'%s\n\n',['% Automagically generated on ',datestr(now)]);
            for ii=1:n
                func_i=func2str(obj(ii).func);
                left_brack=find(func_i=='[',1,'first');
                right_brack=find(func_i==']',1,'last');
                if ~isempty(left_brack)
                    func_i=func_i(left_brack+1:right_brack-1);
                end
                if ~strcmp(func_i(end),';')
                    func_i=strcat(func_i,';'); 
                end
                semicols=strfind(func_i,';');
                [rows,cols]=ind2sub(obj(ii).size,obj(ii).indices);
                siz=num2str(obj(ii).size(1));
                for isiz=2:numel(obj(ii).size)
                    siz=strcat(siz,',',num2str(obj(ii).size(isiz)));
                end
                fprintf(fid,'%s\n',['z=zeros(',siz,');']);
                prev=0;
                for irow=1:numel(rows)
                    next=semicols(irow);
                    rr=int2str(rows(irow));
                    cc=int2str(cols(irow));
                    fprintf(fid,'%s\n',['z(',rr,',',cc,')=',func_i(prev+1:next)]);
                    prev=next;
                end
                if ii==1
                    fprintf(fid,'%s\n\n','Q=z;');
                else
                    fprintf(fid,'%s\n\n','Q=kron(Q,z);');
                end
            end
            fclose(fid);
        end
    end
end