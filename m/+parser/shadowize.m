function [dic,dyn,stat,defs,shadow_tvp,shadow_complementarity]=shadowize(dic,M,eqtn_type)

% define pings
%--------------
dynamic=eqtn_type==1;

definitions=eqtn_type==2;

tvp=eqtn_type==3;

mcp=eqtn_type==4;

planner=eqtn_type==6;

sstate=eqtn_type==5|eqtn_type==7;

dic0=dic;

dic.definitions=struct('name',dic.definitions);

teipes={'endogenous','parameters','exogenous','definitions'};

[ypxd,pat]=preparse(dic,teipes);

dic=dic0;

incid=dic.lead_lag_incidence.before_solve;

% get rid of useless info
%-------------------------
M=[M(:,1),M(:,end)];

M=recollect(M,dynamic);

% Make a copy of the original model
%-----------------------------------
M0=M;

% Generalities + replace switching parameters
%---------------------------------------------
rsp=@(x)parser.add_auxiliary_switching_parameters(x,dic.parameters);

% Non-dynamic equations
M(~dynamic,1)=shadowizer(M(~dynamic,1),ypxd,pat,rsp);

% fast steady state equations
M(dynamic,2)=shadowizer(M(dynamic,2),ypxd,pat,rsp);

% dynamic equations
M(dynamic,1)=shadowizer(M(dynamic,1),ypxd,pat,rsp,incid);

% replace ss with y in all but the dynamic model
%-----------------------------------------------
M(~dynamic,1)=strrep(M(~dynamic,1),'ss','y');

M(dynamic,2)=strrep(M(dynamic,2),'ss','y');

% substitute definitions in all but the definitions
%---------------------------------------------------
shdef=M(definitions,1);

M(~definitions,1)=parser.substitute_definitions(M(~definitions,1),shdef);

M(dynamic,2)=parser.substitute_definitions(M(dynamic,2),shdef);

% now do choppings
%------------------
dyn=do_dynamic();

defs=do_definitions();

do_tvp()

shadow_complementarity=do_mcp();

stat=do_sstate();

do_planner()


    function dyn=do_dynamic()
        
        dyn=struct();
        
        dyn.model=M0(dynamic,1);
        
        dyn.shadow_model=M(dynamic,1);
        
    end


    function defs=do_definitions()
        
        defs=struct();
        
        defs.shadow=M(definitions,1);
        
        defs.original=M0(definitions,1);
        
    end


    function do_tvp()
        
        shadow_tvp=M(tvp,1);
        
    end


    function shadow_complementarity=do_mcp()
        
        shadow_complementarity=M(mcp,1);
        
    end


    function stat=do_sstate()
        
        stat=struct();
        
        stat.model=M0(sstate,1);
        
        stat.shadow_model=M(sstate,1);
        
        stat.shadow_steady_state_model=M(sstate,1);
        
        stat.shadow_fast_ssmodel=M(dynamic,2);
        
    end


    function do_planner()
        
        shadow_plan_syst=M(planner,1);
        
        if isfield(dic.planner_system,'shadow_model')
            
            if ~isempty(shadow_plan_syst)
                
                error('several planner objectives?')
                
            end
            
        else
            
            dic.planner_system.shadow_model=shadow_plan_syst;
            
        end
        
    end


end


function [y,pat]=preparse(dic,teipes)
    
y=struct();

nt=numel(teipes);

list=cell(1,nt);

mytypes=struct('endogenous','y',...
    'parameters','param',...
    'exogenous','x',...
    'definitions','def');

for itype=1:nt
    
    teipe=teipes{itype};  
    
    myteipe=mytypes.(teipe);
    
    sublist={dic.(teipe).name};
    
    ny=numel(sublist);
    
    for ii=1:ny
        
        v=sublist{ii};
        
        y.(v).pos=ii;
        
        y.(v).type=myteipe;
        
    end
    
    list{itype}=sublist(:).';
    
end

list=[list{:}];

pat=parser.cell2matize(list);

end


function M=shadowizer(M,y,pat,rsp,incid)

if nargin<5
    
    incid=[];
    
end

is_dynamic=~isempty(incid);

% replace variables, parameters and definitions in one go in order to avoid
% potential problems with atoms named y,x,param,def...
%------------
rep2=@seek_and_destroy;%#ok<NASGU>

ppat=['\<',pat,'\>'];

if is_dynamic
    
    ppat=[ppat,'(\{[^\}]+\})?'];
    
    M=regexprep(M,ppat,'${rep2($1,$2)}');
    
else
    
    % steady state redux: encode but do not touch the time subscripts
    %--------------------
    M=regexprep(M,ppat,'${rep2($1)}');
    
end

% steady states: do this here in order to avoid any problem with a variable
% being named ss
%-------------------------------------------------------------------------
sspat='\<steady_state\((\w+)\)';

rep1=@re_rep1_; %#ok<NASGU>

M=regexprep(M,sspat,'${rep1($1)}');

% switching parameters to frwz
%-----------------------------
M=rsp(M);

    function o=seek_and_destroy(v,l)
        
        if nargin<2
            
            l=[];
            
        end
        
        loc=y.(v).pos;
        
        t=y.(v).type;
        
        if is_dynamic && strcmp(t,'y')
            
            if isempty(l)
                
                l='0';
                
            else
                
                l=l(2:end-1);
                
            end
            
            c=abs(str2double(l)-2);
            
            loc0=loc;
            
            loc=incid(loc0,c);
            
        end
        
        o=sprintf('%s(%d)',t,loc);
        
    end

    function o=re_rep1_(v)
        
        o=sprintf('ss(%d)',y.pos.(v));
        
    end

end


function M1=recollect(M0,dynamic)

nr=numel(M0);

M1=M0;

for ii=1:nr
    
    Mi=M0{ii};
    
    if isempty(Mi)
        
        continue
        
    end
    
    nc=size(Mi,2);
    
    for jj=1:nc
        
        mij=Mi{2,jj};
        
        if ~isempty(mij) && mij~=0
            
            mij=sprintf('%0.0g',mij);
            
            if mij(1)~='-'
                
                mij=['+',mij]; %#ok<AGROW>
                
            end
            
            Mi{1,jj}=[Mi{1,jj},'{',mij,'}'];
            
        end
        
    end
    
    M1{ii}=cell2mat(Mi(1,:));
    
end

% Take care of fast steady state
%--------------------------------
Md=M1(dynamic,:);

no=cellfun(@isempty,Md(:,2));

Md(no,2)=Md(no,1);

M1(dynamic,:)=Md;

end