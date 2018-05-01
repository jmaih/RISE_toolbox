function string=substitute_definitions(string,definitions)
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
%    See also:

n = nargin;
if n==1
    definitions=string;
end
if ~iscellstr(definitions)
    error('definitions should be a cellstr')
end
ndef=numel(definitions);

defcell=substitute_definitions_in_definitions();

if n==1
    string=defcell(:,2);
else
    for idef=1:ndef
        string=strrep(string,defcell{idef,1},['(',defcell{idef,2},')']);
    end
end

    function defcell=substitute_definitions_in_definitions()
        % get rid of the semicolon
        definitions=strrep(definitions,';','');
        defcell=cell(ndef,2);
        for ii=1:ndef
            def_=definitions{ii};
            equality=strfind(def_,'=');
            defcell(ii,:)={def_(1:equality-1),def_(equality+1:end)};
        end
        for ii=1:ndef
            for jj=ii+1:ndef
                defcell{jj,2}=strrep(defcell{jj,2},defcell{ii,1},['(',defcell{ii,2},')']);
            end
        end
    end
end
