function RawFile=read_file(FileName,definitions)
if nargin<2
    definitions=struct();
end
% definitions=struct('hiring',true,'DMP',true);
% this would be the equivalent of dynare's macroprocessor
RawFile=cell(0,1);
if isempty(FileName)
    return
end
SPACE_DELIMITERS=char([9:13,32]);
fid = fopen(FileName);
iter=0;
% inclusions=false;
allowed=1;
for_looping=0;
if_control=0;
controls={};
tank=cell(0,3);

while 1
    rawline = fgetl(fid);
    if ~ischar(rawline), break, end
    rawline=remove_comments(rawline);
    iter=iter+1;
    if all(isspace(rawline))
        continue
    end
    
    [tokk,rest_]=strtok(rawline,SPACE_DELIMITERS);
    if strncmp(tokk,'@#',2)
        if length(tokk)==2
            oldtokk=tokk;
            [tokk,rest_]=strtok(rest_,SPACE_DELIMITERS); %#ok<STTOK>
            tokk=[tokk,oldtokk];
        end
    end
    if allowed
        switch tokk
            case '@#if'
                allowed=check_condition(rest_,FileName,definitions);
                if_control=if_control+1;
                controls=[controls,'@#if'];
                continue
            case '@#for'
                for_looping=for_looping+1;
                if for_looping>1
                    error('nested for loops not implemented at this point')
                end
                controls=[controls,'@#for'];
                [index,rest_]=strtok(rest_,SPACE_DELIMITERS); %#ok<STTOK>
                index=['@{',index,'}'];
                in=strfind(rest_,'in');
                rest_=rest_(in+2:end);
                colon=strfind(rest_,':');
                start=eval(rest_(1:colon-1));
                finish=eval(rest_(colon+1:end));
                continue
            case '@#else'
                % now we change the behavior
                allowed=~allowed;
                continue
            case '@#elseif'
                allowed=check_condition(rest_,FileName,definitions);
                continue
            case {'@#end','@#endif','@#endfor'}
                rawline={};
                switch controls{end}
                    case '@#for'
                        for_looping=for_looping-1;
                        [rawline,tank]=assign_tank(tank,index,start,finish);
                    case '@#if'
                        if_control=if_control-1;
                    otherwise
                        error('control ended but never started')
                end
                controls(end)=[];
            case '@#include'
                rawline=include_file(rest_,FileName);
            otherwise
                if for_looping
                    tank=build_tank(tank,rawline,FileName,iter);
                    continue
                else
                    rawline={rawline,FileName,iter};
                end
        end
        RawFile=[RawFile;rawline];
    else
        switch tokk
            case '@#for'
                % entering a new control: make things worse
                allowed=allowed-1;
                for_looping=for_looping+1;
                controls=[controls,'@#for'];
            case '@#else'
                allowed=~allowed;
            case '@#elseif'
                allowed=check_condition(rest_,FileName,definitions);
            case '@#if'
                % entering a new control: make things worse
                allowed=allowed-1;
                if_control=if_control+1;
                controls=[controls,'@#if']; %#ok<*AGROW>
            case {'@#end','@#endif','@#endfor'}
                % exiting a control: improve things
                allowed=allowed+1;
                switch controls{end}
                    case '@#for'
                        for_looping=for_looping-1;
                    case '@#if'
                        if_control=if_control-1;
                    otherwise
                        error('control ended but never started')
                end
                controls(end)=[];
            otherwise
                % don't do anything
        end
    end
    
end

fclose(fid);

if for_looping||if_control
    error('for or if flow control not closed')
end

% if inclusions
%     expanded=['expanded_',FileName(1:end-4),FileName(end-3:end)];
%     try %#ok<TRYNC> % this would fail under parallel computing
%         if exist(expanded,'file')
%             delete(expanded)
%         end
%         fid=fopen(expanded,'w+');
%         for irow=1:size(RawFile,1)
%             fprintf(fid,'%s \n',RawFile{irow,1});
%         end
%         fclose(fid);
%     end
% end

end

function tank=build_tank(tank,rawline,FileName,iter)
tank=[tank;{rawline,FileName,iter}];
end
function rf=include_file(rest_,FileName)
quotes=strfind(rest_,'"');
if isempty(rest_)
    error([mfilename,':: missing quotes after statement @#include in ',FileName,' at line',int2str(iter)])
end
if numel(quotes)~=2
    error([mfilename,':: expecting a pair of quotes in ',FileName,' at line ',int2str(iter)])
end
newfile=strtrim(rest_(quotes(1)+1:quotes(2)-1));
if ~exist(newfile,'file')
    error([mfilename,':: file ',newfile,' not found::  ',FileName,' at line ',int2str(iter)])
end
rf=read_file(newfile);
end

function flag=check_condition(rest,FileName,defs)
fields=fieldnames(defs);
for ifield=1:numel(fields)
    ff=fields{ifield};
    tmp=defs.(ff); %#ok<NASGU>
    eval([ff,'=tmp;'])
end
flag=logical(eval(rest));
if isnan(flag)
    error(['expression ',rest,...
        ' in ',FileName,' at line ',int2str(iter),...
        ' is not a valid logical statement'])
end
end

function [rawline,tank]=assign_tank(tank,index,start,finish)
rawline=cell(0,3);
for ii=start:finish
    tank0=tank;
    tank0(:,1)=strrep(tank0(:,1),index,int2str(ii));
    rawline=[rawline;tank0];
end
tank=cell(0,3);
end

function rawline_=remove_comments(rawline_)
% locate comments
loc_=strfind(rawline_,'%');
if ~isempty(loc_)
    rawline_=rawline_(1:loc_(1)-1);
end
loc_=strfind(rawline_,'//');
if ~isempty(loc_)
    rawline_=rawline_(1:loc_(1)-1);
end
end

%================================
% function RawFile=read_file(FileName)
% % this would be the equivalent of dynare's macroprocessor
% RawFile=cell(0,1);
% if isempty(FileName)
%     return
% end
% SPACE_DELIMITERS=char([9:13,32]);
% tank_open=false;
% tank=cell(0,3);
% fid = fopen(FileName);
% iter=0;
% inclusions=false;
% while 1
%     rawline = fgetl(fid);
%     if ~ischar(rawline), break, end
%     rawline=remove_comments(rawline);
%     iter=iter+1;
%     if all(isspace(rawline))
%         continue
%     end
% 
%     % check whether there are files included by parsing rawline
%     [tokk,rest_]=strtok(rawline,SPACE_DELIMITERS);
%     if tank_open
%         if strcmp(tokk,'@#endfor')
%             rawline=cell(0,3);
%             for ii=start:finish
%                 tank0=tank;
%                 tank0(:,1)=strrep(tank0(:,1),index,int2str(ii));
%                 rawline=[rawline;tank0];
%             end
%             tank_open=false;
%             tank=cell(0,3);
%         else
%             tank=[tank;{rawline,FileName,iter}]; %#ok<*AGROW>
%             continue
%         end
%     elseif strcmp(tokk,'@#include')
%         quotes=strfind(rest_,'"');
%         if isempty(rest_)
%             error([mfilename,':: missing quotes after statement @#include in ',FileName,' at line',int2str(iter)])
%         end
%         if numel(quotes)~=2
%             error([mfilename,':: expecting a pair of quotes in ',FileName,' at line ',int2str(iter)])
%         end
%         newfile=strtrim(rest_(quotes(1)+1:quotes(2)-1));
%         if ~exist(newfile,'file')
%             error([mfilename,':: file ',newfile,' not found::  ',FileName,' at line ',int2str(iter)])
%         end
%         rawline=read_file(newfile);
%         inclusions=true;
%     elseif strcmp(tokk,'@#for')
%         inclusions=true;
%         [index,rest_]=strtok(rest_,SPACE_DELIMITERS);
%         index=['@{',index,'}'];
%         in=strfind(rest_,'in');
%         rest_=rest_(in+2:end);
%         colon=strfind(rest_,':');
%         start=eval(rest_(1:colon-1));
%         finish=eval(rest_(colon+1:end));
%         tank_open=true;
%         continue
%     elseif strcmp(tokk,'@#endfor')
%         error([mfilename,':: preprocessor failed ',FileName,' at line ',int2str(iter)])
%     else
%         rawline={rawline,FileName,iter};
%     end
%     RawFile=[RawFile;rawline];
% end
% fclose(fid);
% 
% if inclusions
%     expanded=['expanded_',FileName(1:end-4),FileName(end-3:end)];
%     try %#ok<TRYNC> % this would fail under parallel computing
%         if exist(expanded,'file')
%             delete(expanded)
%         end
%         fid=fopen(expanded,'w+');
%         for ii=1:size(RawFile,1)
%             fprintf(fid,'%s \n',RawFile{ii,1});
%         end
%         fclose(fid);
%     end
% end
% end
% 
% function rawline_=remove_comments(rawline_)
% % locate comments
% loc_=strfind(rawline_,'%');
% if ~isempty(loc_)
%     rawline_=rawline_(1:loc_(1)-1);
% end
% loc_=strfind(rawline_,'//');
% if ~isempty(loc_)
%     rawline_=rawline_(1:loc_(1)-1);
% end
% end
% 
