function [handle,lineloc,mycall]=update_line(item,handle,mycall)
lineloc=nan;
indx=mycall.prefix_list{end};
if any(isletter(item))||strncmp(handle,indx,4) % do not waste time writing constant
    loc=find(strcmp(item,mycall.fid(:,2)));
    if isempty(loc)
        mycall.fid=[mycall.fid;{handle,item}];
        mycall.lineCount=mycall.lineCount+1;
        lineloc=mycall.lineCount;
    else
        handle=mycall.fid{loc,1};
        lineloc=loc;
    end
else
    test=eval(item);
    if ~isnan(test); % instead put the constant to the handle
        handle=sprintf('%0.10g',test);
    end
end
end
