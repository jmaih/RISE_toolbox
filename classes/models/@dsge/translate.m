function outList=translate(obj,inList,order)
% translate -- Translates RISE codes into comprehensible atoms
%
% ::
%
%   outList=translate(obj,inList)
%   outList=translate(obj,inList,order)
%
% Args:
%
%    obj (rise | dsge): scalar model object.
%
%    inList (char | cellstr): List of atoms to translate e.g. y_3,
%       param(4), ss_20, def_10
%
%    order (integer|0|{1}): If order>1 the returned list of atoms is a list
%       of kroneckers in which the number of elements in each items is the
%       order of the kronecker. This is useful for instance to understand
%       what combination of variables makes up a column in a
%       differentiation. If order=0, the translation is with respect to the
%       static, rather than the dynamic model
%
% Returns:
%    :
%
%    - **outList** [cellstr]: List of translated atoms
%
% Example:
%    :
%
%    list=obj.routines.symbolic.probs_times_dynamic{2}
%    outList=translate(obj,list)
%    outList=translate(obj,list,3)
%

if isempty(obj)
    
        outList=cell(0,4);
        
    return
    
end

if nargin<3
    
    order=1;
    
end

if ischar(inList)
    
    inList=cellstr(inList);
    
end

xList=get(obj,'exo_list');

pList=get(obj,'par_list');

yList=get(obj,'endo_list');

dList=get(obj,'def_list');% <---regexprep(obj.definitions.dynamic,'=.+','');

outList=inList;

incidence=obj.lead_lag_incidence.after_solve;%(ov,:)

expr=parser.cell2matize(parser.input_list);

expr=['(?<atom>',expr,')_?\(?(?<digit>\d+)\)?'];

p=regexp(inList,expr,'names');

p=[p{:}];

n=numel(p);

for ii=1:n
    
    t=p(ii).atom;
    
    v=str2double(p(ii).digit);
    
    switch t
        
        case 'y'
            
            outList{ii}=look_up_endogenous();
            
        case 'x'
            
            outList{ii}=look_up_exogenous();
            
        case 'param'
            
            outList{ii}=look_up_parameters();
            
        case 'def'
            
            outList{ii}=look_up_definitions();
            
        case 'ss'
            
            outList{ii}=look_up_sstate();
            
        case 's0'
            
            outList{ii}='first regime';
            
        case 's1'
            
            outList{ii}='second regime';
            
        otherwise
            
            error('Unknown atom')
            
    end
    
end

if order<=1
    
    return
    
end

prev=outList;

for io=2:order
    
    curr=cell(1,n^io);
    
    iter=0;
    
    for ii=1:numel(prev)
        
        for jj=1:n
        
        iter=iter+1;
        
        curr{iter}=[prev{ii},',',outList{jj}];
        
        end
        
    end
    
    prev=curr;
    
end

outList=curr;

    function name=look_up_sstate()
        
        name=[yList{v},'{sstate}'];
        
    end

    function name=look_up_definitions()
        
        name=dList{v};
        
    end

    function name=look_up_parameters()
        
        name=pList{v};
        
    end

    function name=look_up_exogenous()
        
        name=xList{v};
        
    end

    function name=look_up_endogenous()
        % Leads, current, lags
        % v
        
        if order==0
            
            name=yList{v};
            
            return
            
        end
        
        [yr,yc]=find(incidence==v);
        
        name=yList{yr};
        
        if yc==1
            
            name=[name,'{+1}'];
            
        elseif  yc==3
            
            name=[name,'{-1}'];
            
        end
            
    end
    
end