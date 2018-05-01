function [listing,nlist]=insert_definitions(listing,flag)
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

nlist=size(listing,1);

if flag
    isdef=false(nlist,1);
    def_unfinished=false;
    for irow=1:nlist
        thiseq=listing{irow,2};
        thiseq(isspace(thiseq))=[];
        %--------------
        isdef(irow)=strcmp(thiseq(1),'#')||def_unfinished;
        if isdef(irow)            
            if def_unfinished
                rhs=[rhs,thiseq]; %#ok<AGROW>
            else
                eq_sign=strfind(thiseq,'=');
                lhs=thiseq(2:eq_sign-1);
                rhs=thiseq(eq_sign+1:end);
            end
            if strcmp(rhs(end),';')
                def_unfinished=false;
                repstr=['(',rhs(1:end-1),')'];
                expr=['(?<!\w+)',lhs,'(?!\w+)'];
                listing(irow+1:end,2)=regexprep(listing(irow+1:end,2),expr,repstr);
            else
                def_unfinished=true;
            end
        end
        %--------------
    end
    listing=listing(~isdef,:);
    nlist=size(listing,1);
end