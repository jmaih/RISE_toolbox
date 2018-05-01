function xout=code2vector(xcell,devect)
% code2vector - transforms a set of function handles into a single function
%
% ::
%
%
%   xout=code2vector(xcell)
%
% Args:
%
%    - **xcell** [fhandle|cell|struct]:
%      - if cell, it is assumed that all entries are anonymous functions
%        sharing the same inputs. One anonymous function is returned in a cell
%      - if fhandle, the same input is returned
%      - if struct
%          - derivatives are identified by the fields
%            'size','functions','map','partitions'
%          - eval items are identified by the fields 'code','argins','argouts'
%          - else it is assumed we have transition matrix, in which case the
%              output is the same as the input
%
%    - **devect** [true|{false}]: de-vectorize functions.
%
% Returns:
%    :
%
%    - **xout** : vectorized function handle
%
% Note:
%
%    - The routine checks whether the input has ':' or not to writes a
%      consistent unique function
%
% Example:
%
%    See also: UTILS.CODE.CODE2FUNC


xout=xcell;
if isempty(xcell)
    return
end

if nargin<2
    devect=false;
end

if iscell(xcell)
    [xout]=main_engine();
elseif isstruct(xcell)
    derivative_fields={'size','functions','map','partitions'};
    eval_fields={'code','argins','argouts'};
    if all(isfield(xcell,derivative_fields))
        [xout]=derivative_engine();
    elseif all(isfield(xcell,eval_fields))
        [xout]=eval_engine();
    else
        % then it must be a transition matrix
        [xout]=transition_matrix_engine();
    end
elseif isa(xcell,'function_handle')
    xout=xcell;
else
    error('first input must be a cell, a structure or a function handle')
end


    function [xout]=transition_matrix_engine()
        % things are already vectorized... more or less
        xout=xcell;
    end

    function tmp=derivative_engine()
        tmp=xcell;
        order=numel(tmp);
        xout={};
        for io=1:order
            do_one_order(io);
        end
        tmp=rmfield(tmp,'map'); % vectorizer will be used at evaluation...
        
        function do_one_order(oo)
            xcell=tmp(oo).functions;
            tmp(oo).functions=main_engine();
        end
    end

    function [xout]=main_engine()
        n=numel(xcell);
        xout=xcell;
        entry_gate='';
        has_colon=false;
        for item=1:n
            if ~isempty(xcell{item})
                if ~isa(xcell{item},'function_handle')
                    error('all elements in xcell should be function handles')
                end
                xout{item}=func2str(xcell{item});
                if devect
                    xout{item}=devectorize(xout{item});
                end
                if isempty(entry_gate)
                    right_parenth=find(xout{item}==')',1,'first');
                    entry_gate=xout{item}(1:right_parenth);
                end
                if ~has_colon
                    has_colon=any(xout{item}==':');
                end
            end
        end
        xout=strrep(xout,entry_gate,'');
        if has_colon
            xout=strcat(xout(:).',';');
            tp='';
        else
            xout=strcat(xout(:).',',');
            tp='.''';
        end
        xout=cell2mat(xout);
        xout=['[',xout(1:end-1),']',tp];
        xout={str2func([entry_gate,xout])};
    end

    function [code]=eval_engine()
        % eval is already vectorized...
        code=xcell;
    end

    function code=devectorize(code)
        code=strrep(code,',:','');
        code=strrep(code,'./','/');
        code=strrep(code,'.^','^');
        code=strrep(code,'.*','*');
    end
end
