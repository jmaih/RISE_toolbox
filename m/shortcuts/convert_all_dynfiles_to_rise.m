%  convert_all_dynfiles_to_rise Converts Multiple Dynare Files to RISE Format
%  
%  This function searches for Dynare files (.dyn and .mod) in the current
%  directory and its subdirectories. It then converts these Dynare files
%  into RISE files (.rs).  
%  
%  Syntax :
%  
%    convert_all_dynfiles_to_rise()
%  
%    convert_all_dynfiles_to_rise(skip_done) 
%  
%  Inputs :
%  
%    - skip_done [empty|true|{false}]: Optional flag to skip conversion for
%      folders where RISE files already exist.  
%  
%  Outputs :
%  
%    - failures [cell array]: A cell array of failures (if any) that
%      occurred during the conversion. 
%  
%  **Usage**:
%  
%  To convert all Dynare files in the current directory and its
%  subdirectories to RISE format, use: 
%  
%  .. code-block:: matlab
%  
%    convert_all_dynfiles_to_rise()
%  
%  To skip conversion for folders where RISE files already exist, use:
%  
%  .. code-block:: matlab
%  
%    convert_all_dynfiles_to_rise(true)
%  
%  **Details**:
%  
%  - This function searches for Dynare files in the current directory and
%    its subdirectories. 
%          - It then attempts to convert each Dynare file into RISE format.
%          - If a conversion fails for any reason, the function records the
%            failure and continues with the next file. 
%          - The function returns a cell array of failures (if any).
%              - If the `skip_done` flag is set to `true`, the function will
%                skip folders where RISE files already exist. 
%              - Converted RISE files are placed in a folder named
%                'rise_version' in the current directory. 
%  
%              **Example**:
%  
%              .. code-block:: matlab
%  
%              % Convert all Dynare files to RISE format and skip folders
%              with existing RISE files 
%  
%              convert_all_dynfiles_to_rise(true)
%