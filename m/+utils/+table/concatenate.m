function B=concatenate(data,precision)
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

if nargin<2
    precision=[];
end
if isempty(precision)
    precision='%8.6f';
else
    if ~ischar(precision)
        error('precision must be of the type %8.6f')
    end
    precision(isspace(precision))=[];
    if ~strncmp(precision(1),'%',1)||...
            ~all(isstrprop(precision([2,4]),'digit'))||...
            ~isstrprop(precision(end),'alpha')
        error('precision must be of the type %8.6f')
    end
end

B=cell(2,0);
[span,nargs]=size(data);
for icol=1:nargs
    A=parser.any2str(data{1,icol});
    for jrow=2:span
        add_on=data{jrow,icol};
        if iscellstr(add_on)
            add_on=char(add_on);
        end
        if isempty(add_on)
            A=char(A,'--');
        elseif isnumeric(add_on)
            A=char(A,num2str(add_on,precision));
        elseif ischar(A)
            A=char(A,add_on);
        else
            error([mfilename,':: unknown type'])
        end
    end
    B=[B,{A,size(A,2)+2}']; %#ok<AGROW>
end

end
% function B=concatenate(data,precision)
%

% for ii=1:nargs
%     A='';
%     for jrow=1:span
%         add_on=data{jrow,ii};
%         if isnumeric(add_on)
%             add_on=num2str(add_on,precision);
%         end
%         A=char(A,add_on);
%     end
%     A=A(2:end,:);
%     B=[B,{A,size(A,2)+2}']; %#ok<AGROW>
% end