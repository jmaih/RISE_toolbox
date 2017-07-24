function dynaressfile2rise(dynFileName,param_names,riseFileName)

if nargin<3
    
    riseFileName=[];
    
end

if isempty(riseFileName)
    
    riseFileName=[regexprep(dynFileName,'(\w+)\.?\w+?','$1_rise'),'.m'];
    
end

p='p';

newp='newp';

% read input file
%-----------------
raw_code = read_file(dynFileName);

param_names_=cell2mat(strcat(param_names,'|'));

% make a list of the parameters updated in the steady state file
%----------------------------------------------------------------
newpList = parameters_solved_in_ssfile();
% prepare a version to be listed in the RISE file
newpList_=cell2mat(strcat('''',newpList,''','));

% replace parameters computed in steady state
%---------------------------------------------
raw_code=replace_params(raw_code,newpList,newp);

% replace remaining parameters
%------------------------------
oldpList=param_names-newpList;
raw_code=replace_params(raw_code,oldpList,p);

endOfLine=char(10);% sprintf('\n')

% split the code
raw_code=regexp(raw_code,endOfLine,'split');

riseFileName_=regexprep(riseFileName,'\..+','');

header = sprintf('%% Conversion of Dynare file [%s] into RISE file [%s]\n',...
    dynFileName,riseFileName);

% Create and save RISE code.
%---------------------------
timestamp = datestr(now);

header = [header,sprintf('\n %% Done %s.',timestamp)];

raw_code=[
    {['function [y,newp,retcode]=',riseFileName_,'(obj,y,',p,',d,id)']
    ''
    ['% ',riseFileName_,' --  computes the steady state of ... analytically']
    '%'
    '% Syntax'
    '% -------'
    '% ::'
    '%'
    ['%   [y,newp,retcode]=',riseFileName_,'(obj,y,p,d,id)']
    '%'
    '% Inputs'
    '% -------'
    '%'
    '% - **obj** [rise|dsge]: model object (not always needed)'
    '%'
    '% - **y** [vector]: endo_nbr x 1 vector of initial steady state'
    '%'
    '% - **p** [struct]: parameter structure'
    '%'
    '% - **d** [struct]: definitions'
    '%'
    '% - **id** [vector]: location of the variables to calculate'
    '%'
    '% Outputs'
    '% --------'
    '%'
    '% - **y** []: endo_nbr x 1 vector of updated steady state'
    '%'
    '% - **newp** [struct]: structure containing updated parameters if any'
    '%'
    '% - **retcode** [0|number]: return 0 if there are no problems, else return'
    '%   any number different from 0'
    '%'
    '% More About'
    '% ------------'
    '%'
    '% - this is new approach has three main advantages relative to the previous'
    '%   one:'
    '%   - The file is valid whether we have many regimes or not'
    '%   - The user does not need to know what regime is being computed'
    '%   - It is in sync with the steady state model'
    '%'
    '% Examples'
    '% ---------'
    '%'
    '% See also:'
    ''
    header
    ''
    'retcode=0;'
    ''
    'if nargin==1'
    '% list of endogenous variables to be calculated'
    '%----------------------------------------------'
    'y=get(obj,''endo_list'');'
    '% list of parameters to be computed during steady state calculation'
    '%-------------------------------------------------------------------'
    ['newp={',newpList_(1:end-1),'};']
    'else'
    ''
    'newp=struct();'
    ''
    ''}
    raw_code(:)
    {''}
    {'y(id)=ys;'}
    {''}
    {'end'}
    {''}
    {'end'}
    ];

parser.write2file(raw_code,riseFileName);

    function newcode=replace_params(oldcode,whichParams,prefix)
        
        whichParams=cell2mat(strcat(whichParams(:)','|'));
        
        % replace parameters computed in steady state
        
        newcode = regexprep(oldcode,['\<(',whichParams(1:end-1),')\>'],[prefix,'.$1']);
        
    end

    function newpList=parameters_solved_in_ssfile()
        
        param_names_=cell2mat(strcat(param_names,'|'));
        
        % make a list of the parameters updated in the steady state file
        %----------------------------------------------------------------
        newpList = regexp(raw_code,['\<(',param_names_(1:end-1),')\>(\s*=)'],'match');
        
        newpList=cellfun(@(x)x(~isspace(x)),newpList,'uniformOutput',false);
        
        newpList=strrep(newpList,'=','');
        
        newpList=unique(newpList);

    end

end

function raw_code = read_file(dynFileName)

fid = fopen(dynFileName,'r');

if fid == -1
    
    if exist(dynFileName,'file') == false
        
        error('Unable to find ''%s''.',dynFileName);
        
    else
        
        error('Unable to open ''%s'' for reading.',dynFileName);
        
    end
    
end

raw_code = transpose(fread(fid,'*char'));

fclose(fid);

end


