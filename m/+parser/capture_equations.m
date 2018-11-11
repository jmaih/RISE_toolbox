function [block,dic]=capture_equations(dic,listing,block_name)
% INTERNAL FUNCTION
%

if isempty(listing)
    
    block=cell(0,5);
    
    return
    
end

delimiters=parser.delimiters();

y={dic.endogenous.name};

x={dic.exogenous.name};

p={dic.parameters.name};

d=dic.definitions;

yxpd=parser.cell2matize([y(:).',x(:).',p(:).',d(:).']);

tmpdt='t*[+-]*\d+|t';

listing(:,2)=precleaning(listing(:,2),yxpd,tmpdt);

nrows=size(listing,1);

iter=0;

straight=@(x)x(~isspace(x));

eqtn=initialize_equation();

eqtns=eqtn(1:0);

while iter<nrows
    
    iter=iter+1;
    
    iline_=listing{iter,1}; rawline_=straight(listing{iter,2}); file_name_=listing{iter,3};
    
    [model,rest]=strtok(rawline_,delimiters);
    
    if strcmp(model,block_name)
        
        rawline_=rest;
        
    end
    
    if isempty(rawline_)
        
        continue
        
    end
    
    semicol=strfind(rawline_,';');
    
    rest='';
    
    if ~isempty(semicol)
        
        rest=rawline_(semicol+1:end);
        
        rawline_=rawline_(1:semicol);
        
    end
    
    if isempty(eqtn.eqtn)
        
        eqtn=initialize_equation(iline_,file_name_);
        
    end
    
    store_eqtn()
    
end

% parameters in use
%------------------
dic=set_used_atoms(dic,eqtns);

% matlab functions
%-----------------
% matlab_functions(dic)

[eqtns,dic]=replace_definitions(eqtns,dic);

% after the definitions have been inserted, now is the time to set leads
% and lags, not before!!!
eqtns=set_leads_lags(eqtns);

block=eqtns2table(eqtns,yxpd,dic);

dic.exogenous_list={dic.exogenous.name};

dic.parameters_list={dic.parameters.name};

    function store_eqtn()
        
        eqtn.eqtn=[eqtn.eqtn,rawline_];
        
        if isempty(semicol)
            
            return
            
        end
        
        eqtn.finish_line=iline_;
        
        [eqtn,dic,yxpd]=set_type(eqtn,dic,yxpd);

        eqtn=validate_eqtn(eqtn,yxpd,tmpdt,block_name);
        
        eqtn=flush_left_right(eqtn,block_name);
        
        if isempty(eqtns)
            
            eqtns=eqtn;
            
        else
            
            eqtns(end+1)=eqtn;
            
        end
        
        if ~isempty(rest)
            
            eqtn=initialize_equation(iline_,file_name_);
            
            rawline_=rest;
            
            store_eqtn()
            
        else
            
            eqtn=initialize_equation();
            
        end
        
    end

    function eqtns=set_leads_lags(eqtns)
        
        for ieqtn=1:numel(eqtns)
            
            eqtns(ieqtn)=set_one(eqtns(ieqtn));
            
        end
        
        function eqtn=set_one(eqtn)
            
            ll=regexp(eqtn.eqtn,['\<(?<atom>',yxpd,')(\{)(?<val>',tmpdt,')(\})'],'names');
            
            llss=regexp(eqtn.sseqtn,['\<(?<atom>',yxpd,')(\{)(?<val>',tmpdt,')(\})'],'names');
            
            ll=[ll,llss];
            
            if isempty(ll)
                
                return
                
            end
            
            vals=strrep({ll.val},'t','');
            
            zz=cellfun(@isempty,vals);
            
            vals(zz)={'0'}; % arises from expressions such as R(t) or R{t}
            
            vals=cellfun(@str2double,vals,'uniformOutput',true);
            
            eqtn.max_lag=min(0,min(vals));
            
            eqtn.max_lead=max(0,max(vals));
            
            for ii=1:numel(ll)
                
                vn=ll(ii).atom;
                
                [t,vloc]=identify_type(vn);
                
                vv=vals(ii);
                
                if vv<0
                    
                    dic.(t)(vloc).max_lag=min(dic.(t)(vloc).max_lag,vv);
                    
                elseif vv>0
                    
                    dic.(t)(vloc).max_lead=max(dic.(t)(vloc).max_lead,vv);
                    
                end
                
            end
            
            function [t,loc]=identify_type(str)
                
                loc=strcmp(str,y);
                
                if any(loc)
                    
                    t='endogenous';
                    
                else
                    
                    loc=strcmp(str,x);
                    
                    if any(loc)
                        
                        t='exogenous';
                        
                    else
                        
                        loc=strcmp(str,x);
                        
                        t='parameters';
                        
                    end
                    
                end
                
            end
            
        end
        
    end

    function eqtn=initialize_equation(sline,fn)
        
        if nargin==0
            
            sline=nan;
            
            fn='';
            
        end
        
        eqtn=struct('max_lag',0,'max_lead',0,...
            'eqtn','',...
            'sseqtn','',...
            'type','normal',...
            'filename',fn,...
            'is_def',false,... % definitions
            'is_tvp',false,... % endogenous probabilities
            'is_mcp',false,... % complementarity condition
            'start_line',sline,...
            'finish_line',nan);
        
    end

end

function eqtn=flush_left_right(eqtn,block_name)

if ~(eqtn.is_def || eqtn.is_tvp || eqtn.is_mcp || ...
        strcmp(block_name,'steady_state_model')) 
    
    eqtn.eqtn=flusher(eqtn.eqtn);
    
    eqtn.sseqtn=flusher(eqtn.sseqtn);
    
end

    function e=flusher(e)
        
        ekl=strfind(e,'=');
        
        if isempty(ekl)
            
            return
            
        end
        % remove semicolon: some equations (e.g. sstate form) don't have it
        % anyway but make sure to put it back...
        e=strrep(e,';','');
        
        lhs=e(1:ekl-1);
        
        rhs=e(ekl+1:end);
        
        if length(lhs)>length(rhs)
            
            e=[lhs,'-(',rhs,');'];
            
        else
            
            e=['-(',lhs,')','+',rhs,';'];
            
        end
        
    end

end

function [eqtns,dic]=replace_definitions(eqtns,dic)

if ~dic.definitions_inserted
    
    return
    
end

def_eqtns=eqtns([eqtns.is_def]);

dm=definitions_map();

if isempty(dm)
    
    return
    
end

eqtns=eqtns(~[eqtns.is_def]);

n=numel(eqtns);

alleqts=[{eqtns.eqtn},{eqtns.sseqtn}];

for ii=1:size(dm,1)
        
    alleqts=regexprep(alleqts,['\<',dm{ii,1},'\>'],dm{ii,2});
    
end

[eqtns.eqtn]=deal(alleqts{1:n});

[eqtns.sseqtn]=deal(alleqts{n+1:end});

    function dm=definitions_map()
        
        ndef=numel(def_eqtns);
        
        dm=cell(0,2);
        
        if isfield(dic,'definitions_map')
            
            dm=dic.definitions_map;
                        
        end
        
        ndef0=size(dm,1);
        
        dm=[dm;
            cell(ndef,2)];
        
        for jj=ndef:-1:1
            
            stud=def_eqtns(jj).eqtn;
            
            ekl=strfind(stud,'=');
            
            lhs=stud(1:ekl-1);
            
            rhs=['(',strrep(stud(ekl+1:end),';',''),')'];
            
            loc=jj+ndef0;
            
            dm(loc,:)={lhs,rhs};
            
            dm(loc+1:end,2)=regexprep(dm(loc+1:end,2),['\<',lhs,'\>'],rhs);
            
        end
        
        for jj=1:ndef0
            
            lhs=dm{jj,1};
            
            rhs=dm{jj,2};
            
            dm(ndef0+1:end,2)=regexprep(dm(ndef0+1:end,2),['\<',lhs,'\>'],rhs);
            
        end
        
        dic.definitions_map=dm;
        
    end

end

function T=eqtns2table(eqtns,yxpd,dic)

isendo=classify_atoms_variables();

n=numel(eqtns);

T=cell(n,5);

for ii=1:n
    
    e=eqtns(ii);
    
    e1=e.eqtn;
    
    e2=e.sseqtn;
    
    T(ii,:)={do_it(e1),e.max_lag,e.max_lead,e.type,do_it(e2)};
    
end

    function out=do_it(in0)
        
        out=in0;
        
        if isempty(in0)
            
            return
            
        end
        
        in1=regexprep(in0,['\<',yxpd,'\>(\{[^\}]+\})?'],'%$1$2%');
        
        in1=regexp(in1,'%','split');
        
        mt=cellfun(@isempty,in1,'uniformOutput',true);
        
        in1=in1(~mt);
        
        in1=[in1;cell(size(in1))];
        
        for jj=1:size(in1,2)
            
            cinj=in1{1,jj};
            
            if strcmp(cinj(end),'}')
                
                op=strfind(cinj,'{');
                
                in1{1,jj}=cinj(1:op-1);
                % remove t in t+1 t-1 expressions
                in1{2,jj}=xxxdbl(cinj(op+1:end-1));
                
            elseif isvarname(cinj) % <--all(isstrprop(cinj,'alphanum') && isstrprop(cinj(1),'alpha')
                                
                if isendo.(cinj)
                    
                    in1{2,jj}=0;
                    
                end
                
            end
            
        end
        
        out=in1;
        
    end

    function w=classify_atoms_variables()
        
        y={dic.endogenous.name};
        
        x={dic.exogenous.name};
        
        p={dic.parameters.name};
        
        d=dic.definitions;
        
        atoms=[
            y(:)
            x(:)
            p(:)
            d(:)
            ];
        
        na=numel(atoms);
        
        v=zeros(na,1); v(1:numel(y))=1;
        
        w=cell2struct(num2cell(v),atoms,1);
        
    end

    function o=xxxdbl(x)
        
        if strcmp(x,'t')
            
            o=0;
            
        else
            
            o=str2double(strrep(x,'t',''));
            
        end
        
    end

end

function matlab_functions(dic)

types={'endogenous','exogenous','parameters',...'definitions','chain_names','time_varying_probabilities'
    };

for ii=1:numel(types)
    
    t=types{ii};
    
    names={dic.(t).name};
    
    n=numel(names);
    
    if n==0
        
        continue
        
    end
    
    is_matlab=false(1,n);
    
    for jj=1:n
        
        v=names{jj};
        
        is_matlab(jj)=exist([v,'.m'],'file');
        
    end
    
    badnames=names(is_matlab);
    
    if ~isempty(badnames)
        
        disp(['gentle warning (No need to worry): the following atoms(',t,...
            ') are also matlab functions'])
        
        disp(badnames)
        
    end
    
end

end

function dic=set_used_atoms(dic,eqtns)

types={'parameters','exogenous'};

for ii=1:numel(types)
    
    t=types{ii};
    
    p={dic.(t).name};
    
    l=get_list(p);
    
    loc=locate_variables(l,p);
    
    [dic.(t)(loc).is_in_use]=deal(true);
    
end

    function l=get_list(p)
        
        l=regexp({eqtns.eqtn}.',['\<',parser.cell2matize(p),'\>'],'match');
        
        l=[l{:}];
        
        l=unique(l);
        
    end

end

function validate_mcp(e)

ee=e.eqtn(2:end);

dblekl=strfind(ee,'==');

if ~isempty(dblekl) %#ok<STREMP>
    
    error(['complementarity conditions cannot contain double equalities in ',...
        e.filename,' at line ',sprintf('%0.0f',e.start_line)])
    
end

ineq=regexp(ee,'(>=|<=|<|>|\<lt\>|\<le\>|\<gt\>|\<ge\>)','match');

if isempty(ineq)
    
    error(['complementarity conditions must contain at least one inequality sign in ',...
        e.filename,' at line ',sprintf('%0.0f',e.start_line)])
    
end

end

function dic=validate_tvp(e,dic)

ekl=strfind(e.eqtn,'=');

tokk=e.eqtn(2:ekl-1);

[istp,isdiagonal,chain_name]=parser.is_transition_probability(tokk);

iline_=e.start_line;

file_name_=e.filename;

if ~istp
    
    error(['string ''',tokk,''' in ',file_name_,' at line ',...
        sprintf('%0.0f',iline_),' is not an appropriate name for a ',...
        'time-varying transition probability'])
    
end

if isdiagonal
    
    error(['"',tokk,'" is a diagonal transition probabilitiy. Only ',...
        'off-diagonal elements are allowed. Check ',file_name_,...
        ' at line ',sprintf('%0.0f',iline_)])
    
end

dic.time_varying_probabilities=[dic.time_varying_probabilities,{tokk}];

ch_names={dic.markov_chains.name};

loc_chain=find(strcmp(chain_name,ch_names));

if isempty(loc_chain)
    
    error(['markov chain "',chain_name,'" has not be declared. In ',...
        file_name_,' at line ',sprintf('%0.0f',iline_)])
    
end

chain_status=dic.markov_chains(loc_chain).is_endogenous;

if isnan(chain_status)
    
    dic.markov_chains(loc_chain).is_endogenous=true;
    
elseif ~isequal(chain_status,true)
    
    error(['markov chain "',chain_name,'" was previously found to be ',...
        'exogenous and now is endogenous. In ',file_name_,' at line ',...
        sprintf('%0.0f',iline_)])
    
end

end

function e=validate_eqtn(e,yxpd,tmpdt,block_name)

e=split_sstate_check_equalities(e);

patt=['\<',yxpd,'\>((\{)(',tmpdt,')(\}))?']; % new form with only {
% patt=['\<',yxpd,'\>((\(|\{)(',tmpdt,')(\)|\}))?']; old form with ( and {

myreplace=@repl_engine; %#ok<NASGU>

validator(e.eqtn)

validator(e.sseqtn)

    function validator(s)
        
        if isempty(s),return,end
        
        s=strrep(s,';','');
        
        if ~e.is_tvp
            
            ekl=strfind(s,'=');
            
            if ~isempty(ekl)
                
                s=[s(1:ekl-1),'-(',s(ekl+1:end),')'];
                
            end
            
        end
        
        s1=regexprep(s,patt,'${myreplace($1,$2,$3,$4)}');
        
        if strcmp(block_name,'planner_objective')
            
            s1=regexprep(s1,'\<(commitment|discount)\>','${myreplace($1)}');
            
        elseif strcmp(block_name,'steady_state_model')
            % replace all unknown atoms by something computable
            s1=regexprep(s1,'\<([a-zA-Z]+\w*)\>(?!\()','${myreplace($1)}');
                        
        end
        
        echo_chamber(s1,e)
                
    end

    function r=repl_engine(~,~,~,~)
        
        r='0.5';%<--r=sprintf('%0.6g',rand);
        
    end

end

function echo_chamber(s1,e)

% we need to evaluate things here otherwise we might create some variables
% that conflict with subfunctions, static workspaces, etc.

try
    
    eval([s1,';'])
    
catch me
    
    error([me.message,' between line ',int2str(e.start_line),...
        ' and line ',int2str(e.finish_line),' in "',e.filename,'"'])
    
end

end

function eqtn=split_sstate_check_equalities(eqtn)

e=eqtn.eqtn;

pound=strfind(e,'#');

if isempty(pound)
    
    check_equalities(e)
    
    return
    
end

if numel(pound)>1
    
    error(['Too many # in between line ',int2str(eqtn.start_line),...
        ' and line ',int2str(eqtn.finish_line),' in',eqtn.filename])
    
elseif eqtn.is_def

    error(['Definitions cannot have sstate form between line ',int2str(eqtn.start_line),...
        ' and line ',int2str(eqtn.finish_line),' in',eqtn.filename])
    
elseif eqtn.is_mcp

    error(['Restrictions cannot have sstate form between line ',int2str(eqtn.start_line),...
        ' and line ',int2str(eqtn.finish_line),' in',eqtn.filename])
    
elseif eqtn.is_tvp

    error(['Time-varying probabilities cannot have sstate form between line ',int2str(eqtn.start_line),...
        ' and line ',int2str(eqtn.finish_line),' in',eqtn.filename])
    
end

eqtn.sseqtn=e(pound+1:end-1);

eqtn.eqtn=[e(1:pound-1),';'];

check_equalities(eqtn.eqtn)

check_equalities(eqtn.sseqtn)

    function check_equalities(e)
        
        ekl=strfind(e,'=');
        
        if numel(ekl)>1
            
            error(['Too many "=" signs between line ',int2str(eqtn.start_line),...
                ' and line ',int2str(eqtn.finish_line),' in',eqtn.filename])
            
        end
                
    end

end

function [eqtn,dic,yxpd]=set_type(eqtn,dic,yxpd)

chop=true;

switch eqtn.eqtn(1)
    
    case '#'
        
        eqtn.is_def=true;
        
        eqtn.type='def';
        
        ekl=strfind(eqtn.eqtn,'=');
        
        tokk=eqtn.eqtn(2:ekl-1);
        
        dic.definitions=[dic.definitions;{tokk}];
        
        yxpd=[yxpd(1:end-1),'|',tokk,')'];
        
    case '!'
        
        eqtn.is_tvp=true;
        
        eqtn.type='tvp';
        
        dic=validate_tvp(eqtn,dic);
        
    case '?'
        
        eqtn.is_mcp=true;
        
        eqtn.type='mcp';
        
        validate_mcp(eqtn)
        
    otherwise
        
        chop=false;
end

if chop
    
    eqtn.eqtn=eqtn.eqtn(2:end);
    
end

% ! or ? cannot have leads or lags # can have but only if inserted

end

function L=precleaning(L,yxpd,tpmdt)

% parentheses to curly braces
patt=['\<',yxpd,'\((',tpmdt,'|ss|stst|sstate)\)'];

L=regexprep(L,patt,'$1{$2}');

patt=['\<',yxpd,'\{(stst|sstate)\}'];

L=regexprep(L,patt,'steady_state($1)');

% no more ...
L=strrep(L,'...','');

end
