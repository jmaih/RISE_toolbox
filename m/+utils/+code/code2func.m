function func=code2func(code,argins)
% writes a cell array of strings into a structure with a field containing
% functions, a field with information about the size of the matrix to be
% computed and a field with the number of non-zero elements
if nargin<2
    argins=parser.input_list;
end
if ischar(argins)
    argins=cellstr(argins);
end
if ischar(code)
    code=cellstr(code);
end

% put code in vectorized form
%-----------------------------
code=parser.vectorized_model(code,argins);

argins=cell2mat(strcat(argins(:)',','));
main_string=['@(',argins(1:end-1),')'];

siz=size(code);

func=cell(siz);
for irow=1:siz(1)
    for icol=1:siz(2)
        func{irow,icol}=str2func([main_string,code{irow,icol}]);
    end
end

end

%{
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
is_type1=islogical(one_line)||isa(one_line,'double');

if is_type1
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
            code_i=code{ii};
            code_i(isspace(code_i))=[];
            func=[func,code_i,separator]; %#ok<AGROW>
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
else
    
    func=struct('code',cell2mat(code(:)'),'argins',...
        {args},'argouts',{one_line});
    
    %function finalOutput=code2func(xcode,inputList,outputList)
    %finalOutput=struct('code',cell2mat(xcode(:)'),'argins',...
    %    {inputList},'argouts',{outputList});
    %end
end

%}
