function obj=set_auxiliary_switching_parameters(obj)
% perturbation_type :
% - 'm' or 'maih'
% - 'mw' or 'maih_waggoner'
% - 'frwz' or {'frwz',part_list} where is a cell array of strings

nobj=numel(obj);

if nobj>1
    
    for ii=1:nobj
        
        obj(ii)=set_auxiliary_switching_parameters(obj(ii));
        
    end
    
    return
    
end

if sum(obj.parameters.is_switching)==0
    
    return
    
end

is_forced=true;

if isempty(obj.options.solve_perturbation_type)
    
    obj.options.solve_perturbation_type='m';
    
end

perturbation_type=obj.options.solve_perturbation_type;

partition_list=[];

if iscell(perturbation_type)
    
    partition_list=perturbation_type{2};
    
    perturbation_type=perturbation_type{1};
    
    if ~iscellstr(partition_list) %#ok<ISCLSTR>
        
        error('partition list should be a cell array of strings')
        
    end
    
end

pList=get(obj,'par_list(auxiliary)');

pset=pList(:);

pset=[pset,pset];

core=parser.perturbation_control_param_names(true);

switch perturbation_type
    
    case {'m','maih'}
        
        do_initialization()
        
    case {'mw','maih_waggoner'}
        
        loc=strncmp(core(:,1),'iota_2',6);
        
        core{loc,2}=1;
        
        do_initialization()
        
    case 'frwz'
        
        loc=strncmp(core(:,1),'iota_1',6);
        
        core{loc,2}=0;
        
        % derivatives evaluated at sigma=0
        %---------------------------------
        loc=strncmp(core(:,1),'sig',3);
        
        core{loc,2}=0;
        
        do_initialization()
        
        do_set_partitions()
        
    otherwise
        
        error('unknown type of perturbation')
        
end

obj=setup_calibration(obj,pset,is_forced);

    function do_set_partitions()
        
        if isempty(partition_list)
            
            return
            
        end
        
        newList=strcat('part_',partition_list(:).');
        
        % the partitioned guys are 1
        part_pos=locate_variables(newList,pList);
        
        pset(part_pos,2)={1};
        
    end


    function do_initialization()
        
        core_pos=locate_variables(core(:,1),pList);
        
        pset(core_pos,:)=core;
        
        % steady states are 0: they will be changed when solving,
        % if the perturbation strategy is frwz
        pset(strncmp(pset(:,1),'ss_',3),2)={0};
        
        % partitions are 0
        pset(strncmp(pset(:,1),'part_',5),2)={0};
        
    end

end