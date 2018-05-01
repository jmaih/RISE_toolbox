function obj=do_names(obj,val,type)
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


if nargin==0||isempty(obj)
    
    obj=struct();
    
    return
    
end

nval=numel(val);

name_tex_names=cell(0,2);

iter=0;

for ival=1:nval
    
    v=val{1};
    
    if strcmp(v(1),'"')
        
        if ~strcmp(v(end),'"')
            
            error(['wrong specification of tex name "',v,'"'])
            
        end
        
        if iter==0
            
            error(['list should start with model names at ',v])
            
        end
        
        if ~isempty(name_tex_names{iter,2})
            
            error(['two consecutive tex names at ',v])
            
        end
        
        v=v(2:end-1);
        
        if ~isempty(v)
            
            name_tex_names{iter,2}=v;
            
        end
        
    else
        
        iter=iter+1;
        
        name_tex_names(iter,:)={v,[]};
        
    end
    
    val(1)=[];
    
end

% back-tracking tex name
%-----------------------
for irow=1:iter
    
    if isempty(name_tex_names{irow,2})
        
        name_tex_names{irow,2}=name_tex_names{irow,1};
        
    end
    
end

name_tex_names=sortrows(name_tex_names,1);

names=name_tex_names(:,1)';

tex_name=name_tex_names(:,2)';

n=numel(names);

tmp=struct();

tmp.name=names(:)';

tmp.number=n;

if any(strcmp(type,{'parameters','endogenous'}))
    
    % do nothing
    
elseif strcmp(type,'exogenous')
    
    tmp.is_observed=false(1,n);
    
    tmp.number=[n,0];
    
elseif strcmp(type,'observables')
    
    tmp.is_endogenous=true(1,n);
    
    for iname=1:n
        
        loc=find(strcmp(names{iname},obj.endogenous.name));
        
        if ~isempty(loc)
            
            tex_name{iname}=obj.endogenous.tex_name{loc};
            
        else
            
            loc=find(strcmp(names{iname},obj.exogenous.name));
            
            if isempty(loc)
                
                error(['variable "',names{iname},'" not recognized either as endogenous or exogenous'])
            
            end
            
            tex_name{iname}=obj.exogenous.tex_name{loc};
            
            obj.exogenous.is_observed(loc)=true;
            
            tmp.is_endogenous(iname)=false;
        
        end
        
    end
    
    state_endo=locate_variables(tmp.name,obj.endogenous.name,true);
    
    state_endo(isnan(state_endo))=0;
    
    state_exo=locate_variables(tmp.name,obj.exogenous.name,true);
    
    state_exo(isnan(state_exo))=0;
    
    state_id=state_endo+state_exo*1i;
    
    tmp.state_id=state_id(:).';
    
    obj.exogenous.number=[sum(~obj.exogenous.is_observed),sum(obj.exogenous.is_observed)];
    
    tmp.number=[sum(tmp.is_endogenous),sum(~tmp.is_endogenous)];

else
    
    error(['unknown type "',type,'"'])

end

tmp.tex_name=tex_name(:)';

obj.(type)=tmp;

end
