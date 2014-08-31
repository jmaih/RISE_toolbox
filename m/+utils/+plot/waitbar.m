function waitbar(task,x,varargin)

persistent waitbar_handle init_clock 

switch lower(task)
    case 'init'
        if ~isstruct(x)
            error('at initialization, second argument must be a structure')
        end
        waitbar_handle=waitbar(0,x.message,...
            'Name',x.name,...
            'CreateCancelBtn','delete(gcbf)',...
            varargin{:});
        init_clock=clock;
    case 'close'
        if ishandle(waitbar_handle)
            delete(waitbar_handle)
        end
    case 'update'
        if ishandle(waitbar_handle)
            [hrs,mins,secs]=utils.miscellaneous.estimated_time_of_arrival(init_clock,x);
            msg=sprintf('%0.4f complete. ETA: %g:%02g:%02g',x,hrs,mins,secs);
            waitbar(x,waitbar_handle,msg);
            drawnow
        else
            % handle was deleted, don't do anything
        end
    otherwise
        error(['unknown task "',task,'"'])
end

end
