function string=char(this)
n=numel(this);
for i0=1:n
    obj=this(i0);
    % char itself is already taken care of
    if isa(obj.name,'double')
        thisString=sprintf('%0.10g',obj.name); % <-- thisString=num2str(obj.name,10);
    elseif isempty(obj.args) % variable
        thisString=obj.name; % this should be a char
    else
        func_name=obj.name;
        if ismember(func_name,{'mldivide','mrdivide','mpower','mtimes'})
            func_name=func_name(2:end);
        end
        args_=reprocess_arguments(obj.args);
        switch func_name
            case {'normpdf','normcdf'}
                thisString=strcat(func_name,'(',args_{1},',',args_{2},',',args_{3},')');
            case {'exp','log','log10','sin','asin','sinh','cos','acos','cosh',...
                    'tan','atan','tanh','abs','sqrt','isreal','sign'}
                thisString=strcat(func_name,'(',args_{1},')');
            case {'min','max','gt','lt','ge','le'}
                thisString=strcat(func_name,'(',args_{1},',',args_{2},')');
            case {'uplus'}
                thisString=args_{1};
            case {'uminus'}
                thisString=strcat('-(',args_{1},')');
            case {'plus'}
                thisString=strcat(args_{1},'+',args_{2});
            case {'minus'}
                thisString=strcat(args_{1},'-(',args_{2},')');
            case {'times'}
                thisString=strcat('(',args_{1},')*(',args_{2},')');
            case {'power'}
                thisString=strcat('(',args_{1},')^(',args_{2},')');
            case {'rdivide'}
                thisString=strcat('(',args_{1},')/(',args_{2},')');
            otherwise
                error([func_name,' is undefined for objects of class ',mfilename])
        end
    end
    if n>1
        if i0==1
            cellmode=iscell(thisString);
            if cellmode
                string=cell(n,numel(thisString));
            else
                string=cell(size(this));
            end
        end
        if cellmode
            string(i0,:)=thisString;
        else
            string{i0}=thisString;
        end
    else
        string=thisString;
    end
end
    function args=reprocess_arguments(args)
        for ii=1:numel(args)
            if isnumeric(args{ii})
                args{ii}=sprintf('%0.10g',args{ii}); % <-- args{ii}=num2str(args{ii},10);
            else
%                 check=is_ready(args{ii});
                args{ii}=char(args{ii});
%                 if ~check
%                     args{ii}=strcat('(',args{ii},')');
%                 end
            end
        end
        function flag=is_ready(arg)
            alien=~ismember(func_name,{'minus','times','power','rdivide','ldivide'});
            flag=alien||...
                (strcmp(func_name,'minus') && ii==1)||...
                ismember(func_name,arg(1).ready);
            % there may be many arguments but they,ll all share the
            %                     % same functor
        end
    end
end

%{
function string=char(this)
n=numel(this);
for i0=1:n
    obj=this(i0);
    % char itself is already taken care of
    if isa(obj.name,'double')
        thisString=sprintf('%0.10g',obj.name); % <-- thisString=num2str(obj.name,10);
    elseif isempty(obj.args) % variable
        thisString=obj.name; % this should be a char
    else
        func_name=obj.name;
        if ismember(func_name,{'mldivide','mrdivide','mpower','mtimes'})
            func_name=func_name(2:end);
        end
        args_=reprocess_arguments(obj.args);
        switch func_name
            case {'normpdf','normcdf'}
                thisString=strcat(func_name,'(',args_{1},',',args_{2},',',args_{3},')');
            case {'exp','log','log10','sin','asin','sinh','cos','acos','cosh',...
                    'tan','atan','tanh','abs','sqrt','isreal','sign'}
                thisString=strcat(func_name,'(',args_{1},')');
            case {'min','max','gt','lt','ge','le'}
                thisString=strcat(func_name,'(',args_{1},',',args_{2},')');
            case {'uplus'}
                thisString=args_{1};
            case {'uminus'}
                if strcmp(args_{1},'0')
                    thisString='0';
                else
                    thisString=strcat('-(',args_{1},')');
                end
            case {'plus'}
                if (strcmp(args_{1},'0') && strcmp(args_{2},'0'))
                    thisString='0';
                elseif strcmp(args_{1},'0')
                    thisString=args_{2};
                elseif strcmp(args_{1},args_{2})
                    thisString=strcat('2*(',args_{1},')');
                elseif strcmp(args_{2},'0')
                    thisString=args_{1};
                else
                    thisString=strcat(args_{1},'+',args_{2});
                end
            case {'minus'}
                if (strcmp(args_{1},'0') && strcmp(args_{2},'0'))||...
                        strcmp(args_{1},args_{2})
                    thisString='0';
                elseif strcmp(args_{1},'0')
                    thisString=strcat('-(',args_{2},')');
                elseif strcmp(args_{2},'0')
                    thisString=args_{1};
                else 
                    thisString=strcat(thisString,'-(',args_{2},')');
                end
            case {'times'}
                if strcmp(args_{1},'0') || strcmp(args_{2},'0')
                    thisString='0';
                elseif strcmp(args_{1},'1') && strcmp(args_{2},'1')
                    thisString='1';
                elseif strcmp(args_{1},'1')
                    thisString=args_{2};
                elseif strcmp(args_{2},'1')
                    thisString=args_{1};
                elseif strcmp(args_{1},args_{2})
                    thisString=['(',args_{1},')^2'];
                else
                    thisString=strcat('(',args_{1},')*(',args_{2},')');
                end
            case {'power'}
                if strcmp(args_{2},'0')
                    thisString='1';
                elseif strcmp(args_{2},'1')||strcmp(args_{1},'1')
                    thisString=args_{1};
                else
                    thisString=strcat('(',args_{1},')^(',args_{2},')');
                end
            case {'rdivide'}
                if strcmp(args_{1},'0')
                    thisString='0';
                elseif strcmp(args_{1},args_{2})
                    thisString='1';
                elseif strcmp(args_{2},'1')
                    thisString=args_{1};
                else
                    thisString=strcat('(',args_{1},')/(',args_{2},')');
                end
            otherwise
                error([func_name,' is undefined for objects of class ',mfilename])
        end
    end
    if n>1
        if i0==1
            cellmode=iscell(thisString);
            if cellmode
                string=cell(n,numel(thisString));
            else
                string=cell(size(this));
            end
        end
        if cellmode
            string(i0,:)=thisString;
        else
            string{i0}=thisString;
        end
    else
        string=thisString;
    end
end
    function args=reprocess_arguments(args)
        for ii=1:numel(args)
            if isnumeric(args{ii})
                args{ii}=sprintf('%0.10g',args{ii}); % <-- args{ii}=num2str(args{ii},10);
            else
%                 check=is_ready(args{ii});
                args{ii}=char(args{ii});
%                 if ~check
%                     args{ii}=strcat('(',args{ii},')');
%                 end
            end
        end
        function flag=is_ready(arg)
            alien=~ismember(func_name,{'minus','times','power','rdivide','ldivide'});
            flag=alien||...
                (strcmp(func_name,'minus') && ii==1)||...
                ismember(func_name,arg(1).ready);
            % there may be many arguments but they,ll all share the
            %                     % same functor
        end
    end
end

%}