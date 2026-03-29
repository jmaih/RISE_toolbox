%  `PLOT_PROBABILITIES`: Plot smoothed regime and state probabilities.
% 
%    plot_probabilities(model)
% 
%    plot_probabilities(model, specs)
% 
%    plot_probabilities(model, specs, stateList)
% 
%    plot_probabilities(model, specs, stateList, regimeList)
% 
%    plot_probabilities(model, specs, stateList, regimeList, shadeInfo)
% 
%    plot_probabilities(model, specs, stateList, regimeList, shadeInfo, description)
% 
%  **Inputs**:
%    - `model` (abstvar|dsge): Model object.
%    - `specs` (empty|1 x 2 vector): Default number of rows and columns in a figure.
%    - `stateList` (empty|char|cellstr): List of states to plot.
%    - `regimeList` (empty|char|cellstr): List of regimes to plot.
%    - `shadeInfo` (empty|function handle): Function to add shades to the graphs.
%      e.g. `shadeInfo=@(fig)shade(start_finish, colorstr, fig)`.
%    - `description` (empty|struct): description of the elements to plot
% 
%  **Notes**:
%  - The function plots smoothed regime and state probabilities based on the
%    provided model object.
%  - The function allows customization of the plot layout and options.
%  - If `shadeInfo` is provided, it can be a function handle to add shades
%    to the graphs.
% 
%  **Example**:
%    ```
%    plot_probabilities(model, [3, 3], {'state1', 'state2'}, {'regime1', 'regime2'}, @shadeInfoFunction);
%    ```
%  See Also : shade
%