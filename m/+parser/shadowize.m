function [dic,dyn,stat,defs,shadow_tvp,shadow_complementarity]=shadowize(dic,M,eqtn_type)

y=preparse(dic,'endogenous');

p=preparse(dic,'parameters');

x=preparse(dic,'exogenous');

dic0=dic;

dic0.definitions=struct('name',dic0.definitions);

d=preparse(dic0,'definitions');

incid=dic.lead_lag_incidence.before_solve;

% get rid of useless info
%-------------------------
M=[M(:,1),M(:,end)];

M=recollect(M);

% define pings
%--------------
dynamic=eqtn_type==1;

definitions=eqtn_type==2;

tvp=eqtn_type==3;

mcp=eqtn_type==4;

planner=eqtn_type==6;

sstate=eqtn_type==5|eqtn_type==7;

% Make a copy of the original model
%-----------------------------------
M0=M;

% Generalities + replace switching parameters
%---------------------------------------------
rsp=@(x)parser.add_auxiliary_switching_parameters(x,dic.parameters);

M=sspxd(M,y,x,p,d,rsp);

M(dynamic,:)=do_normal(M(dynamic,:),y,incid);

M(~dynamic,1)=do_others(M(~dynamic,1),y);

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


function y=preparse(dic,teipe)

y=struct();

y.list={dic.(teipe).name};

y.pat=parser.cell2matize(y.list);

n=numel(y.list);

y.pos=cell2struct(num2cell(1:n),y.list,2);

end


function M=do_others(M,y)

% steady state redux: encode but do not touch the time subscripts
%--------------------
ypat=['\<',y.pat,'\>'];

rep6=@(u)rep1_(u,y,'y'); %#ok<NASGU>

M=regexprep(M,ypat,'${rep6($1)}');

end


function M=do_normal(M,y,incid)

no=cellfun(@isempty,M(:,2));

M(no,2)=M(no,1);

% dynamic endogenous
%--------------------
ypat=['\<',y.pat,'\>(\{[^\}]+\})?'];

rep5=@rep2_; %#ok<NASGU>

M(:,1)=regexprep(M(:,1),ypat,'${rep5($1,$2)}');

% steady state redux: encode but do not touch the time subscripts
%----------------------------------------------------------------
ypat=['\<',y.pat,'\>'];

rep6=@(u)rep1_(u,y,'y'); %#ok<NASGU>

M(:,2)=regexprep(M(:,2),ypat,'${rep6($1)}');


    function o=rep2_(v,l)
        
        if isempty(l)
            
            l='0';
            
        else
            
            l=l(2:end-1);
            
        end
        
        loc=y.pos.(v);
        
        c=abs(str2double(l)-2);
        
        index=incid(loc,c);
        
        o=sprintf('y(%d)',index);
        
    end


end


function M=sspxd(M,y,x,p,d,rsp)

% steady states
%---------------
sspat='\<steady_state\((\w+)\)';

rep1=@(u)rep1_(u,y,'ss'); %#ok<NASGU>

M=regexprep(M,sspat,'${rep1($1)}');

% parameters
%------------
ppat=['\<',p.pat,'\>'];

rep2=@(u)rep1_(u,p,'param');%#ok<NASGU>

M=regexprep(M,ppat,'${rep2($1)}');

% exogenous
%-----------
xpat=['\<',x.pat,'\>'];

rep3=@(u)rep1_(u,x,'x');%#ok<NASGU>

M=regexprep(M,xpat,'${rep3($1)}');

% definitions
%-----------
dpat=['\<',d.pat,'\>'];

rep4=@(u)rep1_(u,d,'def'); %#ok<NASGU>

M=regexprep(M,dpat,'${rep4($1)}');

% switching parameters to frwz
%-----------------------------
M=rsp(M);

end


function o=rep1_(v,w,id)

o=sprintf('%s(%d)',id,w.pos.(v));

end


function M1=recollect(M0)

nr=numel(M0);

M1=cell(size(M0));

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

end