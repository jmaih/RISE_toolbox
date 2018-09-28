function convert_all_dynfiles_to_rise(skip_done)

if nargin==0
    
    skip_done=[];
    
end

if isempty(skip_done)
    
    skip_done=false;
    
end

is_dyn_file=@(x)length(x)>3 &&...
    (strcmp(x(end-3:end),'.mod')||strcmp(x(end-3:end),'.dyn'));

rise_folder_name='rise_version';

root=pwd();

tmp=dir;

tmp=tmp([tmp.isdir]);

dirnames={tmp.name};

bad_dirs=cellfun(@(x)any(x=='.'),dirnames,'uniformOutput',true);

dirnames=dirnames(~bad_dirs);

for idir=1:numel(dirnames)
    
    allpaths=genpath(pwd);
    
    allpaths=regexp(allpaths,'([^;]+);?','match');
    
    allpaths=strrep(allpaths,';','');
    
    for ipath=1:numel(allpaths)
        
        thisPath=allpaths{ipath};
        
        convert_in_one_folder(thisPath)
        
    end
    
end
        
cd(root)


    function convert_in_one_folder(thisPath)
        
        cd(thisPath)
        
        if skip_done && exist(rise_folder_name,'dir')
            
            return
            
        end
        
        W=dir(thisPath); % W=what(allpaths{ipath})
        
        dynareFiles={W.name};
        
        flag=cellfun(@(x)is_dyn_file(x),dynareFiles,'uniformOutput',true);
        
        dynareFiles=dynareFiles(flag);
        
        detail=regexp(dynareFiles,'(?<fname>\w+)\.(?<xtens>\w+)','names');
        
        if iscell(detail)
            
            detail=[detail{:}];
            
        end
        
        nfiles=numel(detail);
        
        for iname=1:nfiles
            
            if iname==1
                
                fprintf(1,'\n\n %s : %0.0f files\n',thisPath,nfiles);
                
            end
            
            thisName=detail(iname).fname;
            
            thisExtens=detail(iname).xtens;
            
            dyn_name=[thisName,'.',thisExtens];
            
            rise_name=[thisName,'.rs'];
            
            tic
            
            fprintf(1,'Now converting %s into %s',dyn_name,rise_name);
            
            param_file_name=dynare2rise([thisPath,filesep,dyn_name],rise_name);
            
            fprintf(1,'... done in %0.8f secs\n',toc);
            
            if ~exist(rise_folder_name,'dir')
                
                mkdir(rise_folder_name)
                
            end
            
            movefile(rise_name,rise_folder_name)
            
            movefile([param_file_name,'.m'],rise_folder_name)
            
        end
        
    end

end

% xtens=@(x)regexprep(x,'\w+\.(\w+)','$1');
% rootName=@(x)regexprep(x,'(\w+)\.\w+','$1');
