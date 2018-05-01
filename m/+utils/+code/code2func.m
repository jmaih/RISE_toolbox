function func=code2func(code,argins,do_vectorize)
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
%    See also: UTILS.CODE.CODE2VECTOR

% writes a cell array of strings into a structure with a field containing
% functions, a field with information about the size of the matrix to be
% computed and a field with the number of non-zero elements
if nargin<3
    
    do_vectorize=[];
    
    if nargin<2
        
        argins=parser.input_list;
        
    end
    
end

if isempty(argins)
    
    argins=parser.input_list;
end

if isempty(do_vectorize)
    
    do_vectorize=true;
    
end

if ischar(argins)
    
    argins=cellstr(argins);
    
end

if ischar(code)
    
    code=cellstr(code);
    
end

% put code in vectorized form
%-----------------------------
if do_vectorize
    
    code=parser.vectorized_model(code,argins);
    
end

% turn into functions
%---------------------
argins=cell2mat(strcat(argins(:)',','));

main_string=['@(',argins(1:end-1),')'];

func=cellfun(@(x)str2func([main_string,x]),code,'uniformOutput',false);

end
