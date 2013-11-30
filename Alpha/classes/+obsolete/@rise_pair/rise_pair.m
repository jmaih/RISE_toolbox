classdef rise_pair<handle
    properties
        iter=0 % actual size
        jter=0 % nonzero elements
        v=cell(100,1) % record of pairings
        class_name
    end
    methods
        function obj=rise_pair(class_name)
            if nargin
                obj.class_name=class_name;
            end
        end
        function update(obj,location)
            obj.iter=obj.iter+1;
            if ~isempty(location)
                obj.jter=obj.jter+1;
                if obj.jter==numel(obj.v)
                    obj.v(end+100)={[]};
                end
                obj.v{obj.jter}=obj.iter+location*1i;
            end
        end
        function vv=conclude(obj)
            vv=struct('ncols',obj.iter,...
                'pairs',cell2mat(obj.v(1:obj.jter)));
        end
    end
end