function [c,mycall]=char(obj,wrt_index,mycall,first_call)

if nargin<4
    first_call=false;
    if nargin<3
        mycall=[];
        if nargin<2
            wrt_index=[];
        end
    end
end
cols=':';
if ~isempty(mycall)
    prefix=mycall.prefix_list{1};
    if ~isempty(wrt_index)
        do_it=true;
        if iscell(wrt_index)
            nn=numel(wrt_index);
            if nn==1
                wrt_index=wrt_index{1};
            else
                do_it=false;
            end
        end
        if do_it
            item=mat2str(wrt_index);
            item(isspace(item))=',';
            [cols,~,mycall]=update_line(item,'indx',mycall);
        end
    end
end
nobj=numel(obj);
if nobj>1
    siz_obj=size(obj);
    funcs=[obj.func];
    if min(siz_obj)==1 && isa(funcs,'double')
        loc=find(funcs==1);
        if isempty(loc)
            % search for 0 instead
            if any(funcs~=0)
                error('these elements should all be equal to 0: please report this problem to junior.maih@gmail.com')
            end
            c='0';
        elseif numel(loc)>1
            error('too many 1s: please report this problem to junior.maih@gmail.com')
        else
            c=['bigi_(',sprintf('%0.10g',loc),',',cols,')']; % <--- c=sprintf('%s%0.10g%s','bigi_(',loc,',indx)');
        end
    else
        c=cell(siz_obj);
        for iobj=1:numel(obj)
            wrti=[];
            if ~isempty(wrt_index)
                wrti=wrt_index{iobj};
            end
            [c{iobj},mycall]=char(obj(iobj),wrti,mycall);
        end
    end
else
    func_=obj.func;
    args_=obj.args;
    if isa(func_,'double')
        c=sprintf('%0.10g',func_);
    elseif isempty(args_)
        c=func_;
    else
        for iarg=1:numel(args_)
            if isa(args_{iarg},'double')
                args_{iarg}=sprintf('%0.10g',args_{iarg});
            else
                [args_{iarg},mycall]=char(args_{iarg},wrt_index,mycall);
            end
        end
        tmp=sadiff.neat(func_,args_{:});
        if ~isempty(mycall)
            isOp=is_operation(tmp);
        end
        if isempty(mycall)||~isOp
            def_hanlde=tmp;
        else
            mycall.operCount=mycall.operCount+1;
            def_hanlde=create_handle(mycall.operCount,prefix);
            [def_hanlde,~,mycall]=update_line(tmp,def_hanlde,mycall);
        end
        c=def_hanlde;
    end
end
if first_call
    if isempty(mycall)
        c=remove_unnecessary_parentheses(c);
% % %         c=regexprep(c,'\.(\^|\/|\*)','$1');
    end
end
    function flag=is_operation(str)
        flag=false;
        iter=0;
        operation_alphabet='+-*/^(';
        while iter<length(str)
            iter=iter+1;
            if any(str(iter)==operation_alphabet)
                flag=true;
                break
            end
        end
    end
    function out=auxiliary_variable(zzz)
        start='xx_';
        closing='_';
        if isa(zzz,'double')
            % create
            out=[start,sprintf('%0.10g',zzz),closing];
        else
            ns=length(start);
            % check
            out=(zzz(end)==closing) && ...
                length(zzz)>=5 && ...
                strncmp(zzz,start,ns) && ...
                all(ismember(zzz(ns+1:end-1),('0':'9')));
        end
    end
    function flag=contains_letters(str)
        flag=false;
        iter=0;
        alphabet=['a':'z','A':'Z'];
        while iter<length(str)
            iter=iter+1;
            if any(str(iter)==alphabet)
                flag=true;
                break
            end
        end
    end
end
