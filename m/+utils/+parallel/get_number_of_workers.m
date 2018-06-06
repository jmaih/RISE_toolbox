function n=get_number_of_workers()
% INTERNAL FUNCTION: get number of labs available for parallel processing
%
% ::
%
%   n=get_number_of_workers()
%
% Args:
%
%    none
%
% Returns:
%    :
%
%    - **n** [integer]: number of parallel workers available
%

v=ver;
Names={v.Name};
v=v(strcmp(Names,'MATLAB'));
v=str2double(v.Version);
n=0;
success=false;
if v>8.1
    try %#ok<TRYNC>
        pool = gcp('nocreate');
        if ~isempty(pool)
            n=pool.NumWorkers;
        end
        success=true;
    end
end
if ~success
    try %#ok<TRYNC>
        n=matlabpool('size'); %#ok<DPOOL>
    end
end
end
