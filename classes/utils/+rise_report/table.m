% TABLE Table.
%    Tables are used to collect heterogeneous data and metadata into a single
%    container.  Tables are suitable for storing column-oriented or tabular data
%    that are often stored as columns in a text file or in a spreadsheet.  Tables
%    can accommodate variables of different types, sizes, units, etc.  They are
%    often used to store experimental data, with rows representing different
%    observations and columns representing different measured variables.
% 
%    Use the TABLE constructor to create a table from variables in the MATLAB
%    workspace.  Use the readtable function to create a table by reading data
%    from a text or spreadsheet file.
% 
%    The TABLE constructor can also be used to create tables without
%    providing workspace variables, by providing the size and variable
%    types.
% 
%    Tables can be subscripted using parentheses much like ordinary numeric
%    arrays, but in addition to numeric and logical indices, you can use a
%    table's variable and row names as indices.  You can access individual
%    variables in a table much like fields in a structure, using dot
%    subscripting.  You can access the contents of one or more variables using
%    brace subscripting.
% 
%    Tables can contain different kinds of variables, including numeric, logical,
%    character, categorical, and cell.  However, a table is a different class
%    than the variables that it contains.  For example, even a table that
%    contains only variables that are double arrays cannot be operated on as if
%    it were itself a double array.  However, using dot subscripting, you can
%    operate on a variable in a table as if it were a workspace variable.  Using
%    brace subscripting, you can operate on one or more variables in a table as
%    if they were in a homogeneous array.
% 
%    A table T has properties that store metadata such as its variable and row
%    names.  Access or assign to a property using P = T.Properties.PropName or
%    T.Properties.PropName = P, where PropName is one of the following:
% 
%    TABLE metadata properties:
%        Description          - A character vector describing the table
%        DimensionNames       - A two-element cell array of character vectors containing names of
%                               the dimensions of the table
%        VariableNames        - A cell array containing names of the variables in the table
%        VariableDescriptions - A cell array of character vectors containing descriptions of the
%                               variables in the table
%        VariableUnits        - A cell array of character vectors containing units for the variables
%                               in table
%        VariableContinuity   - An array containing a matlab.tabular.Continuity value for each table 
%                               variable, specifying whether a variable represents continuous or discrete 
%                               data values. You can assign 'unset', 'continuous', 'step', or 'event' to
%                               elements of VariableContinuity.
%        RowNames             - A cell array of nonempty, distinct character vectors containing names
%                               of the rows in the table
%        UserData             - A variable containing any additional information associated
%                               with the table.  You can assign any value to this property.
% 
%    TABLE methods and functions:
%      Construction and conversion:
%        table              - Create a table from workspace variables.
%        array2table        - Convert homogeneous array to table.
%        cell2table         - Convert cell array to table.
%        struct2table       - Convert structure array to table.
%        table2array        - Convert table to a homogeneous array.
%        table2cell         - Convert table to cell array.
%        table2struct       - Convert table to structure array.
%      Import and export:
%        readtable          - Create a table by reading from a file.
%        writetable         - Write a table to a file.
%        write              - Write a table to a file.
%      Size and shape:
%        istable            - True for tables.
%        size               - Size of a table.
%        width              - Number of variables in a table.
%        height             - Number of rows in a table.
%        ndims              - Number of dimensions of a table.
%        numel              - Number of elements in a table.
%        horzcat            - Horizontal concatenation for tables.
%        vertcat            - Vertical concatenation for tables.
%      Set membership:
%        intersect          - Find rows common to two tables.
%        ismember           - Find rows in one table that occur in another table.
%        setdiff            - Find rows that occur in one table but not in another.
%        setxor             - Find rows that occur in one or the other of two tables, but not both.
%        unique             - Find unique rows in a table.
%        union              - Find rows that occur in either of two tables.
%      Data manipulation and reorganization:
%        summary            - Print summary of a table.
%        addvars            - Insert new variables at a specified location in a table.
%        movevars           - Move table variables to a specified location.
%        removevars         - Delete the specified table variables.
%        splitvars          - Splits multi-column variables into separate variables.
%        mergevars          - Merges multiple variables into one multi-column variable or a nested table.
%        sortrows           - Sort rows of a table.
%        stack              - Stack up data from multiple variables into a single variable.
%        unstack            - Unstack data from a single variable into multiple variables.
%        join               - Merge two tables by matching up rows using key variables.
%        innerjoin          - Inner join between two tables.
%        outerjoin          - Outer join between two tables.
%        rows2vars          - Reorient rows to be variables of output table.
%        inner2outer        - Invert a nested table-in-table hierarchy.
%        ismissing          - Find elements in a table that contain missing values.
%        standardizeMissing - Insert missing data indicators into a table.
%      Computations on tables:
%        varfun             - Apply a function to variables in a table.
%        rowfun             - Apply a function to rows of a table.
% 
%    Examples:
% 
%       % Create a table from individual workspace variables.
%       load patients
%       patients = table(LastName,Gender,Age,Height,Weight,Smoker,Systolic,Diastolic)
% 
%       % Select the rows for patients who smoke, and a subset of the variables.
%       smokers = patients(patients.Smoker == true, {'LastName' 'Gender' 'Systolic' 'Diastolic'})
% 
%       % Convert the two blood pressure variables into a single variable.
%       patients.BloodPressure = [patients.Systolic patients.Diastolic];
%       patients(:,{'Systolic' 'Diastolic'}) = []
% 
%       % Pick out two specific patients by the LastName variable.
%       patients(ismember(patients.LastName,{'Smith' 'Jones'}), :)
% 
%       % Convert the LastName variable into row names.
%       patients.Properties.RowNames = patients.LastName;
%       patients.LastName = []
% 
%       % Use the row names to pick out two specific patients.
%       patients({'Smith' 'Jones'},:)
% 
%       % Add metadata to the table.
%       patients.Properties.Description = 'Simulated patient data';
%       patients.Properties.VariableUnits =  {''  'Yrs'  'In'  'Lbs'  ''  'mm Hg'};
%       patients.Properties.VariableDescriptions{6} = 'Systolic/Diastolic';
%       summary(patients)
% 
%       % Create a new variable in the table from existing variables.
%       patients.BMI = (patients.Weight * 0.453592) ./ (patients.Height * 0.0254).^2
%       patients.Properties.VariableUnits{'BMI'} =  'kg/m^2';
%       patients.Properties.VariableDescriptions{'BMI'} = 'Body Mass Index';
% 
%       % Sort the table based on the new variable.
%       sortrows(patients,'BMI')
% 
%       % Make a scatter plot of two of the table's variables.
%       plot(patients.Height,patients.Weight,'o')
% 
%       % Create tables from text and spreadsheet files
%       patients2 = readtable('patients.dat','ReadRowNames',true)
%       patients3 = readtable('patients.xls','ReadRowNames',true)
% 
%       % Create a table from a numeric matrix
%       load tetmesh.mat
%       t = array2table(X,'VariableNames',{'x' 'y' 'z'});
%       plot3(t.x,t.y,t.z,'.')
% 
%    See also TABLE, CATEGORICAL
%
%    Reference page in Doc Center
%       doc table
%
%    Other functions named table
%
%       codistributed/table    tall/table
%