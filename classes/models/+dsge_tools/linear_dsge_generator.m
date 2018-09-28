function linear_dsge_generator(riseFile,endo_list)
% linear_dsge_generator -- generates a linear constant-parameter dsge model
%
% ::
%
%   linear_dsge_generator(riseFile)
%   linear_dsge_generator(riseFile,endo_list)
%
% Args:
%
%    riseFile (char): name of the rise file to be generated
%
%    endo_list (cellstr | char| numeric|{3}): list or number of endogenous
%      variables 
%
% Returns:
%    :
%
%    nothing
%
% Note:
%
%    - the endogenous variables are denoted by "x" if the list of
%      endogenous variables is not provided
%
%    - the exogenous variables are denoted by "e"
%
%    - the coefficients on leads of endogenous variables by "ap_i_j"
%
%    - the coefficients on current endogenous variables by "a0_i_j"
%
%    - the coefficients on lags of endogenous variables by "am_i_j"
%
%    - the constant terms by "b_i"
%
%    - the coefficients on exogenous variables by "c_i"
%
%    - where "i" denotes the equation, "j" the variable
%

if nargin<2
    
    endo_list=[];
    
end

if isempty(endo_list),endo_list=3; end

if ischar(endo_list),endo_list=cellstr(endo_list); end

if isnumeric(endo_list)
    
    n=endo_list;
    
    endo_list=[];
    
else
    
    n=numel(endo_list);
    
end

mydot=find(riseFile=='.',1,'first');

if isempty(mydot)
    
    riseFile=[riseFile,'.rs'];
    
else
    
    xt=riseFile(mydot+1:end);
    
    if ~ismember(xt,{'dsge','rs','rz'})
        
        error(['wrong file extention "',xt,'"'])
        
    end
    
end

fid = fopen(riseFile,'w+');

endo_list=set_declarations('endogenous',endo_list);

shk_name=set_declarations('exogenous');

iter_time=0;

incr=1000;

psize=incr;

plist=cell(1,psize);

piter=0;

eqtns=cell(n,1);

for ieq=1:n
        
    eqtn=cell(1,3);
    
    set_position('+')
    
    set_position('0')
    
    set_position('-')
    
    iter_time=0;
    
    eqtn=cell2mat(strcat(eqtn,'+'));
    
    constant=sprintf('b_%0.0f',ieq);
    
    set_param(constant)
    
    shk_param=sprintf('c_%0.0f',ieq);
    
    set_param(shk_param)
    
    eqtn=sprintf('%s+%s+%s*%s=0;',eqtn(1:end-1),constant,shk_param,shk_name{ieq});
    
    eqtns{ieq}=eqtn;
    
end

set_declarations('parameters',plist(1:piter))

set_equations()

fclose(fid);

    function set_position(pos)
        
        switch pos
            
            case '+'
                
                s='p';
                
                t='{+1}';
                
            case '0'
                
                s='0';
                
                t='';
                
            case '-'
                
                s='m';
                
                t='{-1}';
                
        end
        
        subpart=cell(1,n);
        
        for v=1:n
            
            vname=endo_list{v};
            
            pname=sprintf('a%s_%0.0f_%0.0f',s,ieq,v);
            
            set_param(pname)
            
            subpart{v}=sprintf('%s*%s%s',pname,vname,t);
            
        end
        
        subpart=cell2mat(strcat(subpart,'+'));
        
        iter_time=iter_time+1;
        
        eqtn{iter_time}=subpart(1:end-1);
        
    end

    function set_param(pname)
        
        piter=piter+1;
        
        if piter>psize
            
            plist=[plist,cell(1,incr)];
            
            psize=psize+incr;
            
        end
        
        plist{piter}=pname;
        
    end

    function outlist=set_declarations(thetype,list)
        
        if nargin<2
            
            list=[];
            
        end
        
        if isempty(list)
            
            list=cell(1,n);
            
            switch thetype
                
                case 'endogenous'
                    
                    name='x';
                    
                case 'exogenous'
                    
                    name='e';
                    
                otherwise
                    
                    error('unrecognized type')
            end
            
            for ie=1:n
                
                list{ie}=sprintf('%s%0.0f',name,ie);
                
            end
            
        end
        
        batch=cell2mat(strcat(list(:).',',@'));
        
        batch=strrep(batch(1:end-2),'@',' ');
        
        fprintf(fid,'%s %s\n\n',thetype,batch);
        
        if nargout
            
            outlist=list;
            
        end
        
    end

    function set_equations()
        
        fprintf(fid,'model\n\n');
        
        for ii=1:n
            
            fprintf(fid,'   %s\n\n',eqtns{ii});
            
        end
        
    end

end