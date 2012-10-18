classdef rise_variable
    properties(SetAccess=protected)
        name
        tex_name
        id
        % value moved under protection because it is necessary to pass it
        % through set_properties in order to update the dependent
        % properties.
        value
        lb
        ub
        missing
% % % % %         type=nan;
    end
    properties
        det_steady_state
        balanced_growth=0;
        risk
        standard_deviation
        minimum
        maximum
        average
    end
    methods
        % set property utility
        function obj=set_properties(obj,varargin)
            nargs=length(varargin);
            if rem(nargs,2)~=0
                error([mfilename,':: arguments should enter by pairs'])
            end
            object_properties=properties(obj);
            Forbidden={'standard_deviation','minimum','maximum','average','missing'};
            for ii=1:nargs/2
                propname=varargin{2*ii-1};
                if isempty(strcmp(propname,object_properties))
                    error([mfilename,...
                        ':: propname is not a valid property of class ',class(obj)])
                end
                if ~strcmp(propname,Forbidden)
                    propval=varargin{2*ii};
                    obj.(propname)=propval;
                    if strcmp(propname,'value')||...
                            strcmp(propname,'lb')||...
                            strcmp(propname,'ub')
                        obj=obj.update_dependents;
                    end
                end
            end
        end
        % constructor
        function obj=rise_variable(name,varargin)
            if nargin
                obj.name=name;
                obj=set_properties(obj,varargin{:});
            end
        end
    end
    methods(Hidden=true)
        function obj=update_dependents(obj)
            if ~isempty(obj.value)
                obj.missing=false(size(obj.value));
                good=~isnan(obj.value);
                obj.missing(~good)=true;
                if any(good)
                    vv=obj.value(good);
                    obj.standard_deviation=std(vv);
                    obj.minimum=min(vv);
                    obj.maximum=max(vv);
                    obj.average=mean(vv);
                end
            end
            % check that the bounds are consistent
           errmsg=check_bounds_consitency(obj.lb,obj.value,obj.ub);
           if ~isempty(errmsg)
               error([mfilename,':: ',errmsg,' for variable ',obj.name])
           end
        end
    end
end

function errmsg=check_bounds_consitency(low,ct,high)
test= isempty(low)+2*isempty(ct)+3*isempty(high);
retcode=0;
errmsg='';
switch test
    case 0
        retcode=consistency_check(high,ct);
        if retcode==0
            retcode=consistency_check(ct,low);
        end
    case 1
        retcode=consistency_check(high,ct);
    case 2
        retcode=consistency_check(high,low);
    case 3
        retcode=consistency_check(ct,low);
end
switch retcode
    case 1
        errmsg='ct < lb or ub<ct or ub<ct';
    case 2
        errmsg='sizes mismatch in CT vs LB or UB vs CT or UB vs CT';
end

    function retcode=consistency_check(big,small)
        retcode=0;
        if isequal(size(big),size(small))
            tmp=big-small;
            good=~isnan(tmp);
            if any(tmp(good)<0)
                retcode=1;%'upper bound < lower bound';
            end
        else
            retcode=2;%'sizes mismatch between lower and upper bounds';
        end
    end
end
% classdef rise_variable
%     properties(SetAccess=protected)
%         name
%         tex_name
%         id
%     end
%     properties
%         value
%     end
%     methods
%         % set property utility
%         obj=set_properties(obj,varargin)
%         % constructor
%         function obj=rise_variable(varargin) 
%             fields={'name','tex_name','id','value'};
%             for ii=1:nargin
%                 item=varargin{ii}; 
%                 if ii==1
%                     item=deblank(item);
%                 end
%                 obj.(fields{ii})=item;
%             end
%             if isempty(obj.tex_name)
%                 obj.tex_name=obj.name;
%             end
%         end
%     end
% end

%classdef rise_variable %< handle
%    % obj=rise_variable(name,id)
%    % obj=rise_variable(name,id,value)
%    properties
%        value
%    end
%    properties(SetAccess=protected)
%        name
%        id
%    end
%    methods
%        function obj=rise_variable(name,id,varargin)
%            obj.name=deblank(name);
%            obj.id=id;
%            if nargin>2
%                obj.value=varargin{1};
%            end
%        end
%    end
%end

