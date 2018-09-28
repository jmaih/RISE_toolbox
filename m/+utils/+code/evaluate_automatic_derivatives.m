function C=evaluate_automatic_derivatives(funcs,order,engine,varargin)
% INTERNAL FUNCTION: evaluate automatic differentiation
%
% ::
%
%   C=evaluate_automatic_derivatives(funcs,order,tall_and_thin,engine,varargin)
%
% Args:
%
%    - **funcs** [cell array]: functions to differentiate
%    - **order** [integer]: order of differentiation
%    - **engine** [@aplanar|@aplanar_]:
%
% Returns:
%    :
%
%    - **C** [cell array]: derivatives
%

% TODO: remove all functional programming functions, not just if_elseif

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

C=engine(desif_elseif_ize(funcs{1}),with_respect_to,inactive,order);

    function objective=desif_elseif_ize(objective)
        
        s0=strcmp(inactive(:,1),'s0');
        
        if ~any(s0)
            % return immediately otherwise the next statement will crash
            return
            
        end
        
        s0=inactive{s0,2};
        
        s1=inactive{strcmp(inactive(:,1),'s1'),2};
        
        for iobj=1:numel(objective)
            
            objective{iobj}=do_one_function(objective{iobj});
            
        end
        
        function objective=do_one_function(objective)
            
            fstr0=func2str(objective);
            
            % remove calls to if_elseif
            %--------------------------
            [start,finish]=regexp(fstr0,'if_elseif','start','end');
            
            if isempty(start),return,end
            
            op=1;
            
            iter=finish+1;
            
            while op
                
                iter=iter+1;
                
                atom=fstr0(iter);
                
                if atom=='('
                    
                    op=op+1;
                    
                elseif atom==')'
                    
                    op=op-1;
                    
                end
                
            end
            
            preamble=fstr0(1:start-1);
            
            ifelsiff_=fstr0(start:iter);
            
            rest=fstr0(iter+1:end);
            
            studs=regexp(ifelsiff_,'s0==\d+&s1==\d+,','split');
            % ignore if_elseif
            studs=studs(2:end);
            % remove the last element in all entries (commas and
            % closing parenthesis)
            studs1=cellfun(@(x)x(1:end-1),studs,'uniformOutput',false);
            
            voila=regexp(ifelsiff_,'s0==\d+&s1==\d+','match');
            
            ping=strcmp(voila,['s0==',int2str(s0),'&s1==',int2str(s1)]);
            
            objective=str2func([preamble,studs1{ping},rest]);
            
        end
    
    end

end