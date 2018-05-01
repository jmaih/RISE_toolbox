function [fresh_options,missing]=setfield(default_options,varargin)
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


nn=length(varargin);
if nn==1 
    if isempty(varargin{1})
        fresh_options=default_options;
        return
    elseif isstruct(varargin{1})
        new_options=varargin{1};
    else
        error([mfilename,':: varargin must be a structure if its length is 1'])
    end
else
    if rem(nn,2)==0
        names=varargin(1:2:end-1);
        values=varargin(2:2:end);
        new_options=cell2struct(values,names,2);
    else
        error([mfilename,':: arguments must come in pairs'])
    end
end
fresh_options=default_options;
missing=[];
if ~isempty(new_options)
    fields=fieldnames(new_options);
    default_fields=fieldnames(fresh_options);
    missing=false(1,numel(fields));
    for ii=1:numel(fields)
        loc=find(strcmpi(fields{ii},default_fields));
        if ~isempty(loc)
            value=new_options.(fields{ii});
%             if ~isempty(value)
                fresh_options.(default_fields{loc})=value;
%             end
        else
            missing(ii)=true;
        end
    end
    missing=fields(missing);
end
