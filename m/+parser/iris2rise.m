function iris2rise(irisFileName,riseFileName,stderr_name)

if nargin<3
    
    stderr_name=[];
    
    if nargin<2
        
        riseFileName=[];
        
    end
    
end


if isempty(riseFileName)
    
    riseFileName=[regexprep(irisFileName{1},'(\w+)\.\w+','$1'),'.rs'];
    
end

if isempty(stderr_name)
    
    stderr_name='std';
    
end

[equivalences,blknames]=get_equivalences();

% read input file
%-----------------
raw_code = read_file();

rise_code=raw_code;

riseFileName = strtrim(riseFileName);

header = sprintf('Conversion of Iris file [%s] into RISE file [%s]\n',...
    irisFileName{1},riseFileName);

% replace !log_variables !allbut with level_variables
%-----------------------------------------------------
pattern='!log_variables\s+!allbut';
rise_code=regexprep(rise_code,pattern,'level_variables');

% do some translation
%---------------------

do_declarations()

for iequiv=1:numel(blknames)
    
    rise_equiv=equivalences.(blknames{iequiv}).rise;
    
    if ~isempty(rise_equiv)
        rise_code = strrep(rise_code,['!',blknames{iequiv}],...
            rise_equiv);
    end
    
end
rise_code = strrep(rise_code,'=#','=');
rise_code = strrep(rise_code,'!!','#');
rise_code = strrep(rise_code,'...','');				   
rise_code=regexprep(rise_code,'(&|\$)(\w+)(?!\$)','\$($2)');

rblks=parser.initialize_blocks();

rblks=rblks(:,1);
        
% replace the remaining exclamation signs... (preparsing)
rise_code = strrep(rise_code,'!','@# ');

% create standard deviation parameters
%-------------------------------------
stdparams=create_std_parameters();

% Create and save RISE code.
%---------------------------
timestamp = datestr(now);

header = [header,sprintf('\n Done %s.',timestamp)];

recreate_code();

    function recreate_code()
        
        rise_code = [
            regexprep(['%',header],'(\n)','$1%'),sprintf('\n\n'),...
            rise_code
            ];
        
        parser.write2file(rise_code,riseFileName);
        
    end


    function raw_code = read_file()
        
        if ischar(irisFileName)
            
            irisFileName={irisFileName};
            
        end
        
        raw_code='';
        
        for ifile=1:numel(irisFileName)
            
            irisFileName{ifile} = strtrim(irisFileName{ifile});
            
            fid = fopen(irisFileName{ifile},'r');
            
            if fid == -1
                
                if exist(irisFileName{ifile},'file') == false
                    
                    error('Unable to find ''%s''.',irisFileName{ifile});
                    
                else
                    
                    error('Unable to open ''%s'' for reading.',irisFileName{ifile});
                    
                end
                
            end
            
            thisCode=transpose(fread(fid,'*char'));
            
            do_substitutions()
            
            if ifile==1
                
                raw_code = thisCode;
                
            else
                
                raw_code = [raw_code,sprintf('\n\n'),...
                    thisCode
                    ]; %#ok<AGROW>
                
            end
            
            fclose(fid);
        end
        
        function do_substitutions()
            
            [posL,posR]=blocks_locator(thisCode,strcat('!',blknames));
            
            loc_subs=find(strcmp(posR,'substitutions'));
            
            if isempty(loc_subs)
                
                return
                
            end
            
            lcode=length(thisCode);
            
            np=numel(posL);
            
            nsubs=numel(loc_subs);
            
            substretches=cell(1,nsubs);
            
            is_good=0;
            
            for ipos=1:np
                
                pp=posL(ipos);
                
                theType=posR{ipos};
                
                if ipos==1
                    % last one
                    stretch=pp:lcode;
                    
                else
                    % first one
                    stretch=pp:oldpp-1;
                    
                end
                
                oldpp=pp;
                
                if ~strcmp(theType,'substitutions')
                    
                    continue
                    
                end
                
                is_good=is_good+1;
                
                substretches{is_good}=thisCode(stretch);
                
                thisCode=[thisCode(1:stretch(1)-1),...
                    thisCode(stretch(end)+1:end)];
                
            end
            
            do_all_substitutions()
            
            function do_all_substitutions()
                
                nospace=@(x)cellfun(@(z)z(~isspace(z)),x,'uniformOutput',false);
                
                for isub=1:nsubs
                    
                    thisStretch=substretches{isub};
                    
                    batch=regexp(thisStretch,...
                        '(?<lhs>\w+)\s*:\s*=(?<rhs>[^;]*);','names');
                    
                    lhs=nospace({batch.lhs});
                    
                    rhs=nospace({batch.rhs});
                    
                    for ibatch=1:numel(lhs)
                        
                        tosubst=['\$\<',lhs{ibatch},'\>\$'];
                        
                        rhs(ibatch+1:end)=regexprep(rhs(ibatch+1:end),...
                            tosubst,['(',rhs{ibatch},')']);
                        
                        thisCode=regexprep(thisCode,...
                            tosubst,['(',rhs{ibatch},')']);
                        
                    end
                    
                end
                
            end
            
        end
        
    end

    function do_declarations()
        
        [posL,posR]=blocks_locator(rise_code,strcat('!',blknames));
        
        lcode=length(rise_code);
        
        np=numel(posL);
        
        for ipos=1:np
            
            pp=posL(ipos);
            
            theType=posR{ipos};
            
            if ipos==1
                % last one
                stretch=pp:lcode;
                
            else
                % first one
                stretch=pp:oldpp-1;
                
            end
            
            oldpp=pp;
            
            if ~equivalences.(theType).is_declaration
                
                continue
                
            end
            
            rise_code=[rise_code(1:stretch(1)-1),...
                do_one_stretch(rise_code(stretch)),...
                rise_code(stretch(end)+1:end)];
            
        end
        
        function batch=do_one_stretch(batch)
            
            % replace 'description' y with y "description"
            %--------------------------------------------
            express_='''([^'']*)''(\s*|,)(\w+)';
            replace='$3 "$1"';
            batch = regexprep(batch,express_,replace);
            
            % replace "...!!..." with "...#..."
            %----------------------------------
            express_='("[^"!]*)(!!)([^"]*")';
            replace='$1#$3';
            batch = regexprep(batch,express_,replace);
            
            % if it is parameters or exogenous, collect the list
        end
        
    end

    function plist=create_std_parameters()
        
        [posL,posR]=blocks_locator(rise_code,rblks);
        
        loc_subs=find(strcmp(posR,'exogenous'));
        
        if isempty(loc_subs)
            
            return
            
        end
        
        lcode=length(rise_code);
        
        np=numel(posL);
        
        nsubs=numel(loc_subs);
        
        substretches=cell(1,nsubs);
        
        is_good=0;
        
        for ipos=1:np
            
            pp=posL(ipos);
            
            theType=posR{ipos};
            
            if ipos==1
                % last one
                stretch=pp:lcode;
                
            else
                % first one
                stretch=pp:oldpp-1;
                
            end
            
            oldpp=pp;
            
            if ~strcmp(theType,'exogenous')
                
                continue
                
            end
            
            is_good=is_good+1;
            
            substretches{is_good}=rise_code(stretch);
            
        end
        
        plist=do_all_stretches();
        
        function plist=do_all_stretches()
            
            substretches=regexprep(substretches,'"[^"]*"','');
            
            substretches=regexprep(substretches,'\<exogenous\>','');
            
            exo_list=regexp(substretches,'\w+','match');
            
            exo_list=[exo_list{:}];
            
            expre=parser.cell2matize(exo_list);
            
            expre=['\<(',expre,')\>'];
            
            repl=[stderr_name,'_$1*$1'];
            
            block_names=cell2mat(strcat(rblks.','|'));
            
            block_names=['\<(',block_names(1:end-1),')\>'];
            
            [start,finish]=regexp(rise_code,block_names,'start','end');
            
            for ii=numel(finish):-1:1
                
                item=rise_code(start(ii):finish(ii));
                
                if ~strcmp(item,'model')
                    
                    continue
                    
                end
                
                pre=rise_code(1:start(ii)-1);
                
                post='';
                
                if ii==numel(finish)
                    
                    middle=rise_code(start(ii):end);
                    
                else
                    
                    middle=rise_code(start(ii):start(ii+1)-1);
                    
                    post=rise_code(start(ii+1):end);
                    
                end
                
                tmp=regexprep(middle,expre,repl);
                
                % remove/comment out potential model tags
                %-----------------------------------------
                tmp=strrep(tmp,'''','%''');
                
                % remove nonlinear equations
                %----------------------------
                tmp=regexprep(tmp,'=\s*#','=');
                
                rise_code=[pre,tmp,post];
                
            end
            
            % create the list of stdev and add them to the model in a new
            % parameter block
            
            plist=strcat(stderr_name,'_',exo_list);
            
            tmp=cell2mat(strcat(plist,'@'));
            
            tmp=tmp(1:end-1);
            
            tmp=strrep(tmp,'@',' ');
            
            rise_code=[rise_code,sprintf('\n\nparameters\n\n'),tmp];
            
        end
        
    end

end

function [posL,posR]=blocks_locator(code,blknames,sort_type)

if nargin<3
    
    sort_type='descend';
    
end

nblk=numel(blknames);

pos=cell(nblk,2);

isfound=false(1,nblk);

for iblk=1:nblk
    
    p=strfind(code,blknames{iblk});
    
    if ~isempty(p)
        
        np=numel(p);
        
        pos{iblk,1}=p(:);
        
        pos{iblk,2}=repmat(strrep(blknames(iblk),'!',''),1,np);
        
        isfound(iblk)=true;
        
    end
    
end

pos=pos(isfound,:);

posL=cell2mat(pos(:,1));

posR=[pos{:,2}];

[posL,tag]=sort(posL,1,sort_type);

posR=posR(tag);

end

function [equivalences,blknames]=get_equivalences()

equivalences={
    'transition_equations','model',false
    'reporting_equations','model',false
    'measurement_equations','model',false
    'equations','model',false
    'transition_shocks','exogenous',true
    'transition_variables','endogenous',true
    'measurement_variables','endogenous',true
    'variables','endogenous',true
    'shocks','exogenous',true
    'measurement_shocks','exogenous',true
    'exogenous_variables','exogenous',true
    'log_variables','log_vars',true
    'parameters','parameters',true
    'dtrends','',false
    'links','',false
    'sstate_update','',false
    'autoexogenise','',false
    'substitutions','',false
    'allbut','',false
    };
tmp=equivalences;

equivalences=struct();

for ieq=1:size(tmp,1)
    
    equivalences.(tmp{ieq,1})=struct('rise',tmp{ieq,2},'is_declaration',tmp{ieq,3});
    
end

blknames=fieldnames(equivalences);
end
