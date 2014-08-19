clear all
clear classes
%% IMPORTANT
% Please remind Junior Maih to include an example connecting
% a model report to a general report
%% Initialization of the report object and the title page
with_title_page=~false;
if with_title_page
    xrep=rise_report.report('name','report_name',...
        'title','This will be the title of the report',...
        'author',{'first author','second author','third author'},...
        'address',{'address 1','address 2','address 3'},...
        'email',{'first@author.no','second@author.no','third@author.no'},...
        'abstract',{'we try to do this','then we try to do that','hope it works'});
else
    xrep=rise_report.report('name','report_name');
end
%% adding a section title
xrep.section('title','title for the brand new section');
%% adding a subsection title
xrep.subsection('title','title for the brand new subsection');
%% adding a subsubsection title
xrep.subsubsection('title','title for the brand new subsubsection');
%% adding a paragraph title
xrep.paragraph('title','title for the brand new paragraph');
%% adding a subparagraph title
xrep.subparagraph('title','title for the brand new subparagraph');
%% adding an enumeration list
xrep.enumerate('items',{'one thing','leads to another','and back'});
%% adding an itemize list
xrep.itemize('items',{'one thing','leads to another','and back'});
%% adding a verbatim block
xrep.verbatim({'you are supposed to render things','as they just are'})
%% adding some text
xrep.text({'This is a text line',' ','and this is another line'})
%% adding a new page
xrep.newpage()
%% adding a page break
xrep.pagebreak()
%% adding a clear page
xrep.clearpage()
%% adding a clear double page
xrep.cleardoublepage()
%% adding a table
table_itself={
    '','col2','col3','col4'
    'row2',rand,rand,rand
    'row3',rand,rand,rand
    'row4',rand,rand,rand
    'row5',rand,rand,rand
    'row6',rand,rand,rand
    'row7',rand,rand,rand
    'row8',rand,rand,rand
    };
xrep.table('title','title of the table','log',table_itself)
%% adding a figure
fighand=figure('name','just say whatever');
x=1:100;
plot(x,log(x),'linewidth',2)
xrep.figure('title','let''s have a figure','name',fighand);
%% publishing the report
publish(xrep)