%  `PLOT_DATA_AGAINST_PROBABILITIES`: Plot data against smoothed state or
%  regime probabilities. 
% 
%    plot_data_against_probabilities(model, type, specs, shadeInfo)
% 
%  **Inputs**:
%    - `model` (abstvar|dsge): Model object.
%    - `type` (empty|{'state'}|char|cell): Type of plot ('state' or
%      'regime') and list of variables to plot. 
%    - `specs` (empty|1 x 2 vector): Default number of rows and columns in a figure.
%    - `shadeInfo` (empty|function handle): Function to add shades to the graphs.
%      e.g. `shadeInfo=@(fig)shade(start_finish, colorstr, fig)`.
%    - `description` (empty|struct): description of the elements to plot
% 
%  **Notes**:
%  - The function plots data against smoothed state or regime probabilities
%    based on the provided model object. 
%  - The function allows customization of the plot layout and options.
%  - If `shadeInfo` is provided, it can be a function handle to add shades
%    to the graphs. 
% 
%  **Example**:
%    ```
%    plot_data_against_probabilities(model, {'state', {'state1', 'state2'}}, [3, 3], @shadeInfoFunction);
%    ```
%