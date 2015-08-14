function n=get_number_of_workers()
% get_number_of_workers -- get number of labs available for parallel
% processing
%
% Syntax
% -------
% ::
%
%   n=get_number_of_workers()
%
% Inputs
% -------
%
% none
%
% Outputs
% --------
%
% - **n** [integer]: number of labs available
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

v=ver;
Names={v.Name};
v=v(strcmp(Names,'MATLAB'));
v=str2double(v.Version);
n=0;
if v>8.1
    pool = gcp('nocreate');
    if ~isempty(pool)
        n=pool.NumWorkers;
    end
else
    n=matlabpool('size'); %#ok<DPOOL>
end
end
