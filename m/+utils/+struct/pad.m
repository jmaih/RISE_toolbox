function c=pad(varargin)
% pad concatenates fields structures horizontally
%
% Syntax
% -------
% ::
%
%   c=pad(c1,c2,...,cn)
%
% Inputs
% -------
%
% - ci : [struct] 
%
% Outputs
% --------
%
% - c : [struct] structure with fields that horizontally concatenate the
%   fields of c1, c2,..., cn
%
% Remarks
% --------
%
% - All structures must have the same fields and the fields must be
%   concatenable.
%
% Examples
% ---------
%
% c1=struct('a',1,'b',3); c2=struct('a',2,'b',3); c3=struct('a',3,'b',3);
% c=utils.struct.pad(c1,c2,c3) 
%
% See also: 

n=nargin;
silent=true;
for ii=1:n
    if ~isa(varargin{ii},'struct')
        error(sprintf('argument # %0.0f must be a structure',ii))
    end
    fields=fieldnames(varargin{ii});
    if ii==1
        c=struct2cell(varargin{ii});
        main_fields=fields;
    else
        pos=locate_variables(fields,main_fields,silent);
        bad=isnan(pos);
        if any(bad)
            disp(fields(bad))
            error(sprintf('the field(s) above belong to structure # %0.0f but are not part of the first structure',ii))
        elseif numel(pos)~=numel(main_fields)
            error('all structures must have the same fields')
        end
        for ipos=1:numel(pos)
            c{pos(ipos)}=[c{pos(ipos)},varargin{ii}.(fields{ipos})];
        end
    end
end

c=cell2struct(c,main_fields,1);


end