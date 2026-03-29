%  RPRT - Class for generating reports in HTML format.
% 
%  Properties:
%    - Title              - Title for the report.
%    - Author             - Author(s) of the report.
%    - Date               - Date of the report (default: current date and time).
%    - documentclass      - Type of document ('article', 'report', 'book', 'proc', 'minimal', 'slides').
%    - orientation        - Page orientation of the report ('portrait' or 'landscape') (default: 'landscape').
%    - points             - Point size of the text ('10pt', '11pt', '12pt') (default: '12pt').
%    - papersize          - Paper size ('a4paper', 'letterpaper', 'a5paper', 'b5paper', 'executivepaper', 'legalpaper') (default: 'letterpaper').
%    - Body               - Body content of the report.
%    - Footnotes          - Placeholder for footnotes for HTML printing.
%    - maketoc            - Flag indicating whether to include a table of contents in the report (default: false).
%    - ChapterNumber      - Counter for chapter numbering.
%    - SectionNumber      - Counter for section numbering.
%    - SubsectionNumber   - Counter for subsection numbering.
%    - SubsubsectionNumber - Counter for subsubsection numbering.
%    - TableNumber        - Counter for table numbering.
%    - EquationNumber     - Counter for equation numbering.
%    - FigureNumber       - Counter for figure numbering.
%    - FootnoteCount      - Counter for footnote numbering.
%    - filesToDelete      - List of files to delete when saving the report.
%    - packages           - Additional LaTeX packages to include in the document.
% 
%  Methods:
%    - rprt               - Constructor for the rprt class.
%    - chapter           - Add a chapter to the report.
%    - section           - Add a section to the report.
%    - subsection        - Add a subsection to the report.
%    - subsubsection     - Add a subsubsection to the report.
%    - paragraph         - Add a paragraph to the report.
%    - subparagraph      - Add a subparagraph to the report.
%    - text              - Add plain text to the report.
%    - enumerate         - Start an enumerated list in the report.
%    - itemize           - Start a bullet point list in the report.
%    - description       - Start a description list in the report.
%    - table             - Add a table to the report.
%    - longtable         - Add a longtable to the report.
%    - figure            - Add a figure to the report.
%    - create_figure     - creates figures from scratch and adds them to the report.
%    - equation          - Add an equation to the report.
%    - verbatim          - Add verbatim text to the report.
%    - quote             - Add a block quote to the report.
%    - quotation         - Add an inline quotation to the report.
%    - footnote          - Add a footnote to the report.
%    - url               - Add a URL or hyperlink to the report.
%    - save              - Save the report to a file.
%    - publish           - Generates an HTML or PDF file for the report.
%    - newpage           - Add a new page to the report.
%    - clearpage         - flushes the report for tables, figures, etc..
%    - cleardoublepage   - flushes the report for tables, figures, etc..
%    - pagebreak         - Add a page break to the report.
%    - plain             - Add plain text without \text{} command to the report.
%    - rise_file         - Add the content of a rise file to the report.
%    - rise_model_parameterization  - Add the parameter values of a dsge/rise model to the report.
% 
%  Example:
%    report = rprt('Sample Report', 'John Doe','orientation','portrait');
%    report.chapter('Introduction');
%    report.section('Section 1');
%    report.paragraph('This is a sample report.');
%    report.figure('plot.png', 'Figure 1: Plot');
%    report.save('report.html');
%
%    Documentation for rprt
%       doc rprt
%
%