function b=text_verbatim_log(obj)

b=obj.log(:);

if isa(obj,'rise_report.verbatim')
    b=['\begin{verbatim}';b(:);'\end{verbatim}'];
end

end