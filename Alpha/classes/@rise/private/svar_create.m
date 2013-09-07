function obj = svar_create(this,varargin)
% This function creates a rise object of the svar type from a rational
% expectations model. this function is called when the user calls svar on a
% rational expectations models.
pass_along_old=true; % pass along the old model with its contents...

if this.is_svar_model
    obj=this;
else
    varobs_names={this.varobs.name};
    % re-dollarize the tex names so that they are not confused with the
    % model variable names.
    varobs_tex_names=strcat('$',{this.varobs.tex_name},'$');
    varobs=[varobs_names;varobs_tex_names];
    varobs=varobs(:)';
    mysvar=struct('model','svar','var',{varobs});
    oldthis=this;
    % push the options in this and collect them afterwards
    this=set_options(this,varargin{:});
    obj=rise(mysvar);
    obj=set_properties(obj,'options');
    if pass_along_old
        obj=[oldthis;obj];
    else
        obj=[this;obj];
    end
end

end

