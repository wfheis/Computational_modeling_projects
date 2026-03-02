% Brendan Bullock, James Reilly, Will Heisler, Tait Kline
% bbullo27_jereil26_wfheis26_tmklin26_cs346_hw4_extension1.m
% CS346 - Homework 4 - Extension 1 from Base Implementation
% Fall, 2025

% This file corresponds to Homework 4's Exercise 2. However, it is an 
% extension of our "floor project". In this new file, there is far more
% flexibility for experimentation. You can contain the infection to begin
% in a certain area of the grid, infection includes healing stages (where
% infectivity decreases as antigens disappear), probabilities of
% infectivity will increase as infection stages (before healing) advance,
% and there is now a chance that cells will die due to their infection. The
% probability of death will be highest at the last stage of "infectious",
% and will decrease as antigens disappear in the "healing" stage. All other
% factors are the same from the base implementation.

% Variable Initialization % 

m = 50; % rows of grid
n = 50; % columsn of grid
x = 1:m; % initialize x as array [1...m] for indexing 
y = 1:n; % initialize y as array [1...n] for indexing 
p_susceptible = 0.8; % probability a cell will be susceptible
p_infectious = 0.2; % probability a cell will be infectious

num_iter = 300; % amount of iterations of simulation loop
neighborhood_size = 4; % Amount of neighbors to check within neighborhood
prob_catch = 0.20; % the base probability to catch the infection

infect_coefficient = .1; % How much increasing antigens amount increase 
                         % infectivity during infection
deathrate = 0.005; % Value to determine the rate at which infected die
are_dead = zeros(m+2,n+2); % An array to hold which cells are dead

% Grid Initialization % 

rand_grid = rand(m, n); % m x n grid of nums between 1 and 0

% bins for probabliities ranging from 0-100%
edges = [0, p_infectious, p_infectious + p_susceptible, 1];
labels = ["infectious","susceptible","immune"]; % labels for each bin

% categorizes the grid based on the bins (and their corresponding labels) 
categories = discretize(rand_grid, edges, 'categorical', labels);

pop_grid = zeros(m, n); % initialize grid used in simulation

% Define state ranges
infected_states = [1,2,3]; % infection states
healing_states = [4,5]; % healing states
immune_states = [6,7,8,9,10]; % immunity states 
max_immune = immune_states(length(immune_states)); % the maximum immune state

% Mask for constrained spawning region
mask = zeros(m, n); 
constrain_rows = 45:50; % rows to constrain
constrain_cols = 1:5; % columns to constrain
mask(constrain_rows,constrain_cols) = 1; % only these cells "active"
%}

% Assign 0 to the population grid if "susceptible"
pop_grid(categories == "susceptible") = 0;
% Assign either values of infected_states to the population grid if "infectious"
pop_grid(categories == "infectious") = ...
    randi([infected_states(1) healing_states(length(healing_states))],...
    sum(categories(:)=="infectious"), 1);
% Assign (equally) nums. values of immune_states if "immune"
pop_grid(categories == "immune") = ...
    randi([immune_states(1) immune_states(length(immune_states))], ...
    sum(categories(:)=="immune"), 1);

pop_grid = pop_grid.*mask; % apply mask to restrict cell spawning

% Add borders to our population grid
ext_grid = zeros(m+2, n+2) + immune_states(1); % initialize extended grid
ext_grid(2:m+1,2:n+1) = pop_grid; % assign population grid within ext. grid

% Simulation Initializations %

% Cell array to track non-extended grids
gridList = cell(1, num_iter);
gridList{1} = pop_grid;

% Cell array to track extended grids
ext_gridList = cell(1, num_iter);
ext_gridList{1} = ext_grid;

% Anonymous func. for neighbor probability of disease spread
neighbor_prob_sum = @(x,y,infected) (infected(x, y+1) + ...
    infected(x+1, y+2) + infected(x+2, y+1) + ...
    infected(x+1, y)) * prob_catch;

% Simulation Loop %

% Gets an array for the infect coefficients during the "healing" phase
max_infect_coefficient = ... % the highest probability for infection
        (infected_states(length(infected_states))-1)*infect_coefficient+1;
healing_coefficients = ... % an array of probabilities of infeciton in healing phase
       linspace(max_infect_coefficient, 1,length(healing_states)+1);
healing_coefficients(1) = 0; % the first index will always be 0

% Gets an array for the death coefficients during the "healing" phase
death_coefficients = linspace(deathrate*... % array for coefficients
    infected_states(length(infected_states)), ...
    deathrate,length(healing_states)+1);
death_coefficients(1) = 0; % the first index will always be 0

% Creates an array in which we can track specific counts of states
infected_counts = zeros(1, num_iter);
dead_counts = zeros(1, num_iter);
sus_counts = zeros(1, num_iter);

% ismember produces an array of one sand zeros so summing gives counts
% infected_counts includes the infecteed state and healing state
infected_counts(1) = sum(ismember(ext_grid, infected_states), "all") ...
    + sum(ismember(ext_grid, healing_states), "all");
% All values that are 0 are susceptible
sus_counts(1) = sum(ismember(ext_grid, 0), "all");
% Sums up all dead which innitialy is 0
dead_counts(1) = sum(are_dead, "all");

ismember(ext_grid, infected_states)

for i = 2:num_iter % for number of iterations,
    pop_grid = gridList{i-1}; % pull previous non-extended grid state
    ext_grid = ext_gridList{i-1}; % pull previous extended grid state
    
    % Gets just the values of cells who are infected
    infected_cells = ismember(ext_grid, infected_states); % logical array
    stages_of_infected = ext_grid.*infected_cells; % removes all non-infected
    
    % Gets just the values of cells who are healing
    healing_cells = ismember(ext_grid, healing_states); % logical array
    stages_of_healing = ext_grid.*healing_cells; % removes all non-healing

    % Calculates the new values for infectivity in "infectious" stage cells, 
    % determined by state and the infect_coefficient, with a min of 1
    infected_coefficient_values = max(0, infected_cells+...
    (stages_of_infected - 1)*infect_coefficient);

    % Calculates the new values for infectivity in "healing" stage cells,
    % determined by their state, with earlier healing stages having
    % infectivity.
    temp_healing_coe = zeros(m+2, n+2)+... 
        (max(1, stages_of_healing-(healing_states(1)-2)));
    final_temp_healing_coe = healing_coefficients(temp_healing_coe);

    % combines both healing and infectious stages coefficients into one
    % array.
    final_infect_coefficient_values = ...
        infected_coefficient_values+final_temp_healing_coe;

    % Calculates the new probabilities of death in "infectious" stage
    % cells.
    death_coefficient_values = max(0, infected_cells+...
    (stages_of_infected)*deathrate-1);

    % Calculates the new probabilities of deaht in "healing" stage cells.
    % Will decrease as antigens disappear and cell heals.
    temp_death_coe = zeros(m+2, n+2)+...
        (max(1, stages_of_healing-(healing_states(1)-2)));
    final_temp_death_coe = death_coefficients(temp_death_coe);
    final_death_coefficient_values = ...
        death_coefficient_values+final_temp_death_coe;
    
    % Determine probabilities of infected grid utilizing anonymous func.
    infect_prob = neighbor_prob_sum(x,y,final_infect_coefficient_values);

    % initialize extended probabliity grid with border padding, assign
    % probabilities.
    infect_prob_pad = zeros(m+2, n+2); 
    infect_prob_pad(2:m+1, 2:n+1) = infect_prob;
    % Determine susceptible cells, assign probabilities they catch disease
    susceptible = (ext_grid == 0); % logical array, suscept. = 1
    prob_to_catch = infect_prob_pad .* susceptible;
    
    % Determine cells that will not be infected in extended and
    % non-extended grids
    safe_from_infection_pad = prob_to_catch == 0; % logical array, 0% = 1
    safe_from_infection = safe_from_infection_pad(2:end-1, 2:end-1); 
    % Determine cells that are exposed to infection in non-extended grid
    exposed_to_infection = ~safe_from_infection; %logical array, >0% = 1
    
    % Apply state advancements % 

    % Add one to all elements, except if they are susceptible to disease
    non_suceptible_state_no_pad = (pop_grid+1) .* (pop_grid ~= 0);
    % Add the first immune state value in the case of border cells
    non_suceptible_state_pad = zeros(m+2, n+2) + immune_states(1);                                                                
    % Assign no border grid to be within the extended grid with borders
    non_suceptible_state_pad(2:m+1, 2:n+1) = non_suceptible_state_no_pad;
    
    rand_grid2 = rand(m+2,n+2); % a random probability between 0 and 1
    % logical array, if probCatch > prob = 1
    will_catch = (prob_to_catch > rand_grid2); 
    % Add one to all susceptible cells that caught the infection
    before_immune_reset = non_suceptible_state_pad + will_catch;

    % Uses the rand_grid2 to check if the cell dies
    will_die = (final_death_coefficient_values > rand_grid2);
    % Updates the are_dead array with new dead cells
    are_dead = are_dead+will_die;

    % If the cell ran out of immunity, reset the cell to susceptible again
    immune_reset = (before_immune_reset > max_immune) == 0;
    % Final grid - all updates applied to all cells
    final = before_immune_reset .* immune_reset;
    % At any corresponding cell in are_dead is not 0, meaning the cell is
    % dead, change its value to -1 in final to indicate the cell is dead
    final(are_dead ~= 0) = -1;
     
    % Update variables %

    ext_grid = final;
    pop_grid = ext_grid(2:m+1, 2:n+1);
    gridList{i} = pop_grid; % update non-extended grid cell array
    ext_gridList{i} = ext_grid; % update extended grid cell array
    
    % Updates the counts for the arrays
    infected_counts(i) = sum(infected_cells + healing_cells, "all");
    dead_counts(i) = sum(are_dead,"all");
    sus_counts(i) = sum(ext_grid == 0,"all");

end

% Visualization of Results %

% A function to visualize infection spread,taking in our gridList,  the 
% interval for change between states, the bounds for the color map, and 
% stage arrays of immunity/infection/healing.
function visualize_infection_spread(gridList, interval,...
    lowerCmapBound, upperCmapBound, infected_states, immune_states, ...
    healing_states)

    % Colormap limits
    upper = upperCmapBound;
    lower = lowerCmapBound;
    % Custom colormap applicable to the simulation
    susceptible_color = [0.5 0.5 0.5];   % 0/susceptible = grey
    dead = [0 0 0];   % -1/dead = black
    
    % immunity color channel arays, with length equal to the number of 
    % immune states.
    immune_cb = linspace(0.4,1,length(immune_states))'; % blue channel
    immune_cg = linspace(0,0.8,length(immune_states))'; % green channel
    immune_cr = linspace(0,0.6,length(immune_states))'; % red channel
    immune_colors = [immune_cr immune_cg immune_cb]; % final immune colors

    % infected color channel arays, with length equal to the number of 
    % infected states.
    infected_cb = linspace(0.6,0.25,length(infected_states))'; % blue channel
    infected_cg = linspace(0.6,0.25,length(infected_states))'; % green channel
    infected_cr = linspace(1,0.8,length(infected_states))'; % red channel
    infected_colors = [infected_cr infected_cg infected_cb]; % final infected colors
    
    % healing color channel arays, with length equal to the number of 
    % healing states. Since infectivity decreases during healing, colors
    % will reflect infectivity values.
    healing_cb = linspace(infected_cb(length(infected_cb)-1)... % blue channel
        ,0.6,length(healing_states))';
    healing_cg = linspace(infected_cg(length(infected_cb)-1)... % green channel
        ,0.6,length(healing_states))';
    healing_cr = linspace(infected_cr(length(infected_cb)-1)... % red channel
        ,1,length(healing_states))';
    healing_colors = [healing_cr healing_cg healing_cb]; % final healing colors
    
    % The custom color map utilizing all defined color cases.
    customMap = [
        dead;
        susceptible_color;
        infected_colors;
        healing_colors;
        immune_colors
    ];
    
    figure; % create new figure
    colormap(customMap); % assign color map

    % Loop through frames and display them
    for i = interval:interval:length(gridList) % for all stored grids
        clf; % clear window
        
        data = gridList{i}; % take frame of interest
        % data % Prints data to confirm the grid is displayed properly
        imagesc(data); % display data with full color map range
        clim([lower upper]); % assign color map range
        axis equal; % set aspect ratio to be equal
        axis tight; % axes fit "tightly" around data
        colorbar; % show colorbar
        
        title(sprintf('Frame: %d', i)); % title 
        fprintf('Press any key for next frame...\n'); % advancement
        waitforbuttonpress; % wait to advance until a button is clicked
    end
end

time = 1:num_iter;

figure;
hold on;
plot(time, infected_counts);
plot(time, dead_counts);
plot(time, sus_counts);

xlabel('Time Step');
ylabel('Number of Cells');
title('Counts of Cell States Over Time');
legend('Infected', 'Dead', 'Susceptible');

% call our function for all of our simulation data.
visualize_infection_spread(gridList, 1, -1, ...
    immune_states(length(immune_states)), infected_states, ...
    immune_states, healing_states)
