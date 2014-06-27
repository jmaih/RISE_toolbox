classdef rise_anonymous
    properties
        size
        indices
        func
        auxiliary=''
        args
        type='fhandle'
    end
    methods
        function obj=rise_anonymous(siz,fhandle,ind,auxiliary)
            if nargin
                obj.size=siz;
                obj.indices=ind;
                if nargin<4
                    auxiliary='';
                end
                if ~isempty(auxiliary)
                    obj.type='eval';
                end
                if isa(fhandle,'function_handle')
                    get_arguments();
                end
                obj.auxiliary=auxiliary;
                obj.func=fhandle;
            end
            function get_arguments()
                string=func2str(fhandle);
                closing_par=find(string==')',1,'first');
                obj.args=regexp(string(3:closing_par-1),',','split');
                if strcmp(obj.type,'eval')
                    fhandle=string(closing_par+1:end);
                else
                    if isempty(strfind(string,'double'))
                        fhandle=[string(1:closing_par),'double(',string(closing_par+1:end),')'];
                        fhandle=str2func(fhandle);
                    end
                end
            end
        end
        function val=eval(obj,varargin)
            n=numel(obj);
            val=cell(1,n);
            for iarg=1:numel(obj(1).args)
                eval([obj(1).args{iarg},'=varargin{iarg};'])
            end
            % it is assumed that all objects in the vector share the same
            % auxiliary equations.
            eval(obj(1).auxiliary)
            for ii=1:n
                val{ii}=zeros(obj(ii).size);
                if strcmp(obj(ii).type,'eval')
                    val{ii}(obj(ii).indices)=eval(obj(ii).func);
                else
                    val{ii}(obj(ii).indices)=obj(ii).func(varargin{:});
                end
            end
        end
        function [Q,retcode]=kron(obj,varargin)
            n=numel(obj);
            Q0=eval(obj,varargin{:}); %#ok<*EVLC>
            Q=Q0{1};
            for ii=2:n
                val_i=Q0{ii};
                if isequal(Q,1)
                    Q=val_i;
                elseif isequal(val_i,1)
                    % no need to update Q
                else
                    Q=kron(Q,val_i);
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
            args_=strcat(obj(1).args,',');
            args_=args_(1:end-1);
            fid=fopen([filename,'.m'],'w');
            fprintf(fid,'%s\n',['function Q=',filename,'(',cell2mat(args_),')']);
            fprintf(fid,'%s\n\n',['% Automagically generated on ',datestr(now)]);
            
            % first write the auxiliary equations
            auxil=obj(1).auxiliary;
            if ~isempty(auxil)
                auxil=strcat(regexp(auxil,';','split'),';');
                auxil=auxil(1:end-1);
                fprintf(fid,'%s\n',' %% Auxiliary equations');
                for iaux=1:numel(auxil)
                    fprintf(fid,'%s\n',auxil{iaux});
                end
            end
            
            fprintf(fid,'%s\n','%% Main equations');
            for ii=1:n
                func_i=obj(ii).func;
                if isa(func_i,'function_handle')
                    func_i=func2str(func_i);
                end
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
                % this does not scale to more than 2 dimensions. I should
                % correct that some time.
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