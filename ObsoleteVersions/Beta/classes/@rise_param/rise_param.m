classdef rise_param
    properties(SetAccess=protected)
        name
        id
        is_switching=false;
    end
    properties
        tex_name
        startval
    end
    methods
        % set property utility
        obj=set_properties(obj,varargin)
        % constructor
        function obj=rise_param(varargin) 
            fields={'name','tex_name','id','startval','is_switching'};%
            nn=nargin;
            if ~(rem(nn,2)==0)
                error([mfilename,':: arguments must come in pairs'])
            end
            for ii=1:nn/2
                ff=varargin{(ii-1)*2+1};
                loc=find(strcmp(ff,fields));
                if isempty(loc)
                    error([mfilename,':: unrecognized field ',ff])
                end
                vv=varargin{(ii-1)*2+2};
                if strcmp(vv,'name')
                    vv=deblank(ff);
                end
                obj.(fields{loc})=vv;
            end
            if isempty(obj.tex_name)
                obj.tex_name=obj.name;
            end
        end
    end
end

%classdef rise_param
%    properties
%        name
%        tex_name
%        id
%        startval
%    end
%    methods
%        % set property utility
%        obj=set_properties(obj,varargin)
%        % constructor
%        function obj=rise_param(name,tex_name,id,startval)
%			for ii=1:nargin
%				obj.(
%			end
%            if nargin<4
%                startval=[];
%                if nargin<2
%                    id=[];
%                    if nargin<1
%                        name='';
%                    end
%                end
%            end
%            obj.name=deblank(name);
%            obj.id=id;
%            obj.startval=startval;
%        end
%    end
%end