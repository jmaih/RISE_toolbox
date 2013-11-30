function func=code2func(code,args,one_line)
% writes a cell array of strings into a structure with a field containing
% functions, a field with information about the size of the matrix to be
% computed and a field with the number of non-zero elements
if ischar(args)
    args=cellstr(args);
end
if ischar(code)
    code=cellstr(code);
end
args=cell2mat(strcat(args(:)',','));
main_string=['@(',args(1:end-1),')'];
siz=size(code);
if one_line
    if min(siz)>1
        error('cannot write this code in one line')
    elseif min(siz)==0
        func=[];
        return
    end
    func=[main_string,'['];
    vertical=siz(2)==1;
    if vertical
        separator=';';
    else
        separator=',';
    end
    for ii=1:max(siz)
        func=[func,code{ii},separator]; %#ok<AGROW>
    end
    func=strrep([func(1:end-1),']'],';]',']');
    func={str2func(func)};
    % set the size to nan so that computation takes place in one call
    % without the need to preallocate a matrix.
    siz=nan;
else
    func=cell(siz);
    for irow=1:siz(1)
        for icol=1:siz(2)
            func{irow,icol}=str2func([main_string,code{irow,icol}]);
        end
    end
end

func=struct('size',{siz},...
    'functions',{func},...
    'nnz_derivs',prod(siz));

