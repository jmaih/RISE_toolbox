function C=evaluate_automatic_derivatives(funcs,order,engine,varargin)
% evaluate_automatic_derivatives -- automatic differentiation
%
% ::
%
%
%   C=evaluate_automatic_derivatives(funcs,order,tall_and_thin,engine,varargin)
%
% Args:
%
%    - **funcs** [cell array]: functions to differentiate
%
%    - **order** [integer]: order of differentiation
%
%    - **engine** [@aplanar|@aplanar_]:
%
% Returns:
%    :
%
%    - **C** [cell array]: derivatives
%
% Note:
%
% Example:
%
%    See also:

input_list=parser.input_list();
ninput=numel(input_list);
if length(input_list)~=ninput
    error('insufficient number of input arguments')
end
prototype_inputs=regexprep(input_list,'(\w+)','$1_');%strcat(input_list,'_');
% remove s0 and s1
prototype_inputs=regexprep(prototype_inputs,'s\d+_','');
prototype_inputs=prototype_inputs(~cellfun(@(x)isempty(x),prototype_inputs));
% replace them by s
prototype_inputs=union(prototype_inputs,'s');
allsymb=parser.collect_symbolic_list(funcs{1},prototype_inputs);
% make sure all the active are part of the symbolic list
with_respect_to=funcs{2};
allsymb=union(allsymb,with_respect_to);

nsymb=numel(allsymb);
symbvals=nan(nsymb,1);
for ii=1:ninput
    inp_name=input_list{ii};
    vlength=length(inp_name);
    pos=strncmp(inp_name,allsymb,vlength);
    if any(pos)
        if any(strcmp(inp_name,{'s0','s1'}))
            symbvals(pos)=varargin{ii};
        else
            locs=regexprep(allsymb(pos),'\w+_(\d+)','$1');
            locs=cellfun(@(x)str2double(x),locs,'uniformOutput',false);
            symbvals(pos)=varargin{ii}(cell2mat(locs));
        end
    end
end
inactive=setdiff(allsymb,with_respect_to);
aloc=locate_variables(with_respect_to,allsymb);
inacloc=locate_variables(inactive,allsymb);
allsymb=[allsymb(:),num2cell(symbvals(:))];
with_respect_to=allsymb(aloc,:);
inactive=allsymb(inacloc,:);

C=engine(funcs{1},with_respect_to,inactive,order);