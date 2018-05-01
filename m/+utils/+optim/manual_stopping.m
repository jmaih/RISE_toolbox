function flag=manual_stopping(action)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if nargin==0
    action=0;
end
flag=0;
ManualStoppingFile='ManualStoppingFile.txt';
switch action
    case 0
        fid=fopen(ManualStoppingFile,'w');
        fclose(fid);
    otherwise
        rawfile = char(textread(ManualStoppingFile,'%s','delimiter','\n','whitespace','','bufsize',40000));
        sizerawfile=size(rawfile);
        if sizerawfile(2)
            iline=0;
            while iline < sizerawfile(1)
                iline=iline+1;
                rawline = rawfile(iline,:);
                rawline(isspace(rawline))=[];
                if ~isempty(rawline)
                    if any(strcmp(rawline,{'progress','report'}))
                        flag=3;
                        utils.optim.manual_stopping(0);
                    elseif any(strcmp(rawline,{'pause','keyboard'}))
                        flag=2;
                    else
                        flag=1;
                    end
                    break
                end
            end
        end
end
end
