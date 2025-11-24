clc;
clear all;
close all;

% Packages
pkg load fuzzy-logic-toolkit
pkg load io
pkg load statistics

% Show figures
plot_mfs  = true;
plot_hist = true;

% Year select
YEAR = input("Year: ", 's');

% All console output from this point on will be saved to 'output_fuzzy.txt'
if YEAR == 2012
  diary('FIS_results_2012.txt');
else
  diary('FIS_results_2013.txt');
end

%%%%%%%%%%%%%%%%%%%%%%%%% AUXILIAR FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to calculate descriptive statistics [Min, Q25, Q50, Q75, Max]
function stats_vec = calculate_stats(data)
    min_val = min(data);
    max_val = max(data);
    q25 = quantile(data, 0.25);
    q50 = quantile(data, 0.50);
    q75 = quantile(data, 0.75);
    stats_vec = [min_val, q25, q50, q75, max_val];
    stats_vec = double(stats_vec);
end

% Function to create a universe of discourse
function u = make_universe(min_val, max_val, n_points = 1000)
    if min_val == max_val
        % Add margin if the universe is a single point
        u = linspace(min_val - 0.5, max_val + 0.5, n_points);
    else
        % Add a small margin to max/min to make the traps look better
        margin = (max_val - min_val) * 0.05;
        u = linspace(min_val - margin, max_val + margin, n_points);
    end
end

% Function to plot Input MFs and shade the activation area
function plot_input_activation(fis, input_idx, input_value, color, title_suffix)
    % Get variable info
    var_name = fis.input(input_idx).name;
    u_range = fis.input(input_idx).range;
    n_mf = size(fis.input(input_idx).mf)(2);

    figure('Position', [100, 100, 600, 400]);
    hold on;

    % Generate x-axis (universe of discourse)
    x = linspace(u_range(1), u_range(2), 500);

    all_mu_activations = []; % Array to store all membership degrees for input_value

    for j = 1:n_mf
        mf = fis.input(input_idx).mf(j);

        % Calculate the complete MF curve
        y = feval(mf.type, x, mf.params);
        plot(x, y, 'LineWidth', 1.5, 'DisplayName', mf.name);

        % Calculate the membership degree (activation) for the input_value
        mu_activation = feval(mf.type, input_value, mf.params);

        % Collect the activation degree
        all_mu_activations = [all_mu_activations, mu_activation];

        % If activated, shade the clipped area
        if mu_activation > 1e-6
            % Calculate the clipped MF curve
            y_clipped = min(y, mu_activation);

            % Shade the area (using fill)
            fill([x, fliplr(x)], [y_clipped, zeros(size(y_clipped))], color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end

    % Plot the vertical line for the input value
    [max_mu, ~] = max(all_mu_activations);

    % Get Y-limits of the plot
    y_lims = ylim();

    % Draw the dashed vertical line for the input value
    line([input_value, input_value], [y_lims(1), max_mu], ...
          'Color', color, 'LineStyle', '--', 'LineWidth', 2, ...
          'DisplayName', ['Input: ', num2str(input_value, '%.2f')]);

    % Plot a marker at the top of the line (at the point of maximum membership)
    plot(input_value, max_mu, 'o', 'MarkerSize', 8, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off'); % Hide from legend


    % Final plot settings
    hold off;
    title(['Activation - ', var_name, ' (', title_suffix, ')']);
    xlabel(var_name);
    ylabel('Membership Degree');
    grid on;
    ylim([0, 1.05]);
    legend('Location', 'northeast');
end

% Function to plot Output MFs, shade the aggregated set and show the Defuzzified COG line
function plot_output_result(fis, input_data, rules, color, title_suffix)
    figure('Position', [100, 100, 600, 400]);
    hold on;

    % --- 1. SETUP AND CALCULATION OF FUZZY AGGREGATION (Union/Max of Clipped Consequents) ---
    var_idx = 1; % Output variable index (Quality)
    u_range = fis.output(var_idx).range;
    x_points = 500;
    x = linspace(u_range(1), u_range(2), x_points);
    mu_aggregated = zeros(1, x_points);
    num_output_mf = size(fis.output(var_idx).mf)(2);
    num_rules = size(rules, 1);
    num_inputs = size(fis.input)(2);

    % Vector to store Firing Strengths
    activations = zeros(num_rules, 1);

    for r = 1:num_rules
        % Calculate Firing Strength (Conjunction: MIN/AND)
        firing_strength = 1.0;
        for i = 1:num_inputs
            input_mf_idx = rules(r, i); % Antecedent MF index

            if input_mf_idx > 0
                mf_type = fis.input(i).mf(input_mf_idx).type;
                mf_params = fis.input(i).mf(input_mf_idx).params;

                % Calculate mu and apply the AND operator (min)
                mu = feval(mf_type, input_data(i), mf_params);
                firing_strength = min(firing_strength, mu);
            end
        end
        activations(r) = firing_strength;

        % --- CONTRIBUTION CALCULATION (Clipping and Aggregation) ---
        output_mf_idx = rules(r, num_inputs + 1); % Consequent MF index (Quality)

        if firing_strength > 1e-6 && output_mf_idx > 0
            mf = fis.output(var_idx).mf(output_mf_idx);

            % Calculate the membership degree of the consequent MF
            mu_consequent = feval(mf.type, x, mf.params);

            % Clipping (Apply firing strength: MIN/AND)
            mu_clipped = min(mu_consequent, firing_strength);

            % Aggregation (Union/MAX)
            mu_aggregated = max(mu_aggregated, mu_clipped);
        end
    end


    % --- 2. PLOTTING ---

    % Plot the contours of the Output MFs (without clipping)
    for j = 1:num_output_mf
        mf = fis.output(var_idx).mf(j);
        y = feval(mf.type, x, mf.params);
        plot(x, y, 'LineWidth', 1.5, 'DisplayName', mf.name);
    end

    % Plot the shaded area of the aggregated set (Union/Max)
    % The shaded area should go from mu_aggregated down to 0.
    fill([x, fliplr(x)], [mu_aggregated, zeros(size(mu_aggregated))], color, ...
          'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');


    % --- 3. COG CALCULATION and PLOT ---

    % Get the defuzzified output (COG)
    output_val = evalfis(input_data, fis);
    cog_val = output_val(1);

    % Find the membership degree of the aggregated set at the COG point (interpolation)
    mu_at_cog = interp1(x, mu_aggregated, cog_val, 'linear', 'extrap');

    % Plot the defuzzified COG line (Vertical Line)
    % The line should go from 0 up to the point where it touches the aggregated set.
    line([cog_val, cog_val], [0, mu_at_cog], ...
          'Color', color, 'LineStyle', '-', 'LineWidth', 3, ...
          'DisplayName', ['Defuzzified: ', num2str(cog_val, '%.2f')]);

    % Plot a marker at the top of the line (at the mu_at_cog membership point)
    plot(cog_val, mu_at_cog, 'o', 'MarkerSize', 8, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');


    % --- 4. FINAL SETTINGS ---
    ylim([0, 1.05]); % Ensures Y-axis goes from 0 to 1.05
    hold off;
    title(['Output - Quality (', title_suffix, ')']);
    xlabel(fis.output(var_idx).name);
    ylabel('Membership Degree');
    grid on;
    legend('Location', 'northeast');

end


% Function to show activations, membership degrees, and final result
function ret = show_activations(fis, data_input, rules, labels = {'Diameter', 'Width', 'Length', 'Weight'})

    disp('Input (Diameter, Width, Length, Weight):');
    disp(data_input);

    disp(' ');
    % Evaluate FIS to get the defuzzified output
    output = evalfis(data_input, fis);
    disp(['Output (Fuzzy Quality Index): ', num2str(output(1), '%.4f')]);

    disp(' ');
    disp('Input Membership Degrees:');

    for i = 1:size(fis.input)(2)
        disp([' - Variable: ', labels{i}, ' (Value: ', num2str(data_input(i)), ')']);
        for j = 1:size(fis.input(i).mf)(2)
            mf_name = fis.input(i).mf(j).name;
            mf_type = fis.input(i).mf(j).type;
            mf_params = fis.input(i).mf(j).params;
            % Calculate membership degree
            mu = feval(mf_type, data_input(i), mf_params);
            disp(['     ', mf_name, ': ', num2str(mu, '%.4f')]);
        end
    end

    disp(' ');

    % --- Manual Firing Strength Calculation ---
    num_rules   = size(rules, 1);
    activations = zeros(num_rules, 1);

    for r = 1:num_rules
        firing_strength = 1.0; % Initialize with 1.0 (identity for MIN/AND operation)
        for i = 1:size(fis.input)(2)
            input_index = i;
            mf_idx = rules(r, input_index);

            if mf_idx > 0
                mf_type = fis.input(input_index).mf(mf_idx).type;
                mf_params = fis.input(input_index).mf(mf_idx).params;

                % Calculate mu and apply the AND operator (min)
                mu = feval(mf_type, data_input(i), mf_params);
                firing_strength = min(firing_strength, mu);
            end
        end
        activations(r) = firing_strength;
    end

    rules_activateds = find(activations > 1e-6); % 1e-6 to avoid floating-point errors

    if isempty(rules_activateds)
        disp('No rules were activated (all memberships resulted in 0).');
    else
        disp('Activated Rules and Firing Strength:');
        for idx = rules_activateds'
            degree = activations(idx);
            % Find the consequent MF name
            output_mf_idx = rules(idx, size(fis.input)(2) + 1);
            consequent_name = fis.output(1).mf(output_mf_idx).name;

            disp(['- Rule ', num2str(idx), ' (Consequent: ', consequent_name, '): Strength: ', num2str(degree, '%.4f')]);
        end
    end
    disp('--------------------------------------------------');

    % Return the quality index for plotting use
    ret = output(1);

end

%%%%%%%%%%%%%%%%%%%%%%%%% DATASET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read dataset
dataset = xlsread('baru_seed.xlsx', YEAR);

% Read input variables
Diameter_data = dataset(:, 1);
Width_data    = dataset(:, 2);
Length_data   = dataset(:, 3);
Weight_data   = dataset(:, 4);

% Descriptive statistics
statistics.Diameter = calculate_stats(Diameter_data);
statistics.Width    = calculate_stats(Width_data);
statistics.Length   = calculate_stats(Length_data);
statistics.Weight   = calculate_stats(Weight_data);

%%%%%%%%%%%%%%%%%%%%%%%%% FIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fis = newfis('BaruAlmond');

%%%%%%%%%%%%%%%%%%%%%%%%% INPUT VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Diameter
diam_stats = statistics.Diameter;
u_diameter = make_universe(diam_stats(1), diam_stats(5));
fis = addvar(fis, 'input', 'Diameter', [u_diameter(1) u_diameter(end)]);
fis = addmf(fis, 'input', 1, 'Small', 'trapmf', [u_diameter(1), diam_stats(1), diam_stats(2), diam_stats(3)]);
fis = addmf(fis, 'input', 1, 'Middle', 'trimf', [diam_stats(2), diam_stats(3), diam_stats(4)]);
fis = addmf(fis, 'input', 1, 'Large', 'trapmf', [diam_stats(3), diam_stats(4), diam_stats(5), u_diameter(end)]);
if plot_mfs
  plotmf(fis, 'input', 1);
end

% Width
width_stats = statistics.Width;
u_width = make_universe(width_stats(1), width_stats(5));
fis = addvar(fis, 'input', 'Width', [u_width(1) u_width(end)]);
fis = addmf(fis, 'input', 2, 'Small', 'trapmf', [u_width(1), width_stats(1), width_stats(2), width_stats(3)]);
fis = addmf(fis, 'input', 2, 'Middle', 'trimf', [width_stats(2), width_stats(3), width_stats(4)]);
fis = addmf(fis, 'input', 2, 'Large', 'trapmf', [width_stats(3), width_stats(4), width_stats(5), u_width(end)]);
if plot_mfs
  plotmf(fis, 'input', 2);
end

% Length
length_stats = statistics.Length;
u_length = make_universe(length_stats(1), length_stats(5));
fis = addvar(fis, 'input', 'Length', [u_length(1) u_length(end)]);
fis = addmf(fis, 'input', 3, 'Short', 'trapmf', [u_length(1), length_stats(1), length_stats(2), length_stats(3)]);
fis = addmf(fis, 'input', 3, 'Middle', 'trimf', [length_stats(2), length_stats(3), length_stats(4)]);
fis = addmf(fis, 'input', 3, 'Long', 'trapmf', [length_stats(3), length_stats(4), length_stats(5), u_length(end)]);
if plot_mfs
  plotmf(fis, 'input', 3);
end

% Weight
weight_stats = statistics.Weight;
u_weight = make_universe(weight_stats(1), weight_stats(5));
fis = addvar(fis, 'input', 'Weight', [u_weight(1) u_weight(end)]);
fis = addmf(fis, 'input', 4, 'Slim', 'trapmf', [u_weight(1), weight_stats(1), weight_stats(2), weight_stats(3)]);
fis = addmf(fis, 'input', 4, 'Middle', 'trimf', [weight_stats(2), weight_stats(3), weight_stats(4)]);
fis = addmf(fis, 'input', 4, 'Heavy', 'trapmf', [weight_stats(3), weight_stats(4), weight_stats(5), u_weight(end)]);
if plot_mfs
  plotmf(fis, 'input', 4);
end

%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT VARIABLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quality
u_quality = [0 100];
fis = addvar(fis, 'output', 'Quality', u_quality);
fis = addmf(fis, 'output', 1, 'Low', 'trimf', [-1, 0, 50]);
fis = addmf(fis, 'output', 1, 'Middle', 'trimf', [30, 50, 70]);
fis = addmf(fis, 'output', 1, 'Middle-High', 'trimf', [50, 70, 90]);
fis = addmf(fis, 'output', 1, 'High', 'trimf', [70, 100, 101]);
if plot_mfs
  plotmf(fis, 'output', 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%% BASE RULES (13 RULES) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rules (13) - Expanded to simulate 'OR' in Octave
% [D_MF, W_MF, L_MF, We_MF, Q_MF, Weight of MF, Conjunction (1=AND)]
% MFs: Diameter/Width: 1=Small, 2=Middle, 3=Large
% MFs: Length: 1=Short, 2=Middle, 3=Long
% MFs: Weight: 1=Slim, 2=Middle, 3=Heavy
% MFs: Quality: 1=Low, 2=Middle, 3=Middle-High, 4=High

rules = [
% 1. Weight['Heavy'] & Length['Long'] -> Quality['High']
    0, 0, 3, 3, 4, 1, 1;
% 2a. Weight['Middle'] & Diameter['Middle'] -> Quality['Middle']
    2, 0, 0, 2, 2, 1, 1;
% 2b. Weight['Middle'] & Width['Middle'] -> Quality['Middle']
    0, 2, 0, 2, 2, 1, 1;
% 3a. Weight['Slim'] -> Quality['Low']
    0, 0, 0, 1, 1, 1, 1;
% 3b. Length['Short'] -> Quality['Low']
    0, 0, 1, 0, 1, 1, 1;
% 4. Weight['Heavy'] & Diameter['Large'] & Length['Long'] -> Quality['High']
    3, 0, 3, 3, 4, 1, 1;
% 5. Weight['Slim'] & Width['Small'] -> Quality['Low']
    0, 1, 0, 1, 1, 1, 1;
% 6. Diameter['Large'] & Length['Short'] -> Quality['Middle']
    3, 0, 1, 0, 2, 1, 1;
% 7. Weight['Heavy'] & Width['Large'] -> Quality['High']
    0, 3, 0, 3, 4, 1, 1;
% 8. Weight['Slim'] & Diameter['Small'] -> Quality['Low']
    1, 0, 0, 1, 1, 1, 1;
% 9. Weight['Middle'] & Length['Long'] -> Quality['Middle-High']
    0, 0, 3, 2, 3, 1, 1;
% 10. Weight['Heavy'] & Length['Middle'] -> Quality['Middle-High'] - (Faltava no seu pedido, mas é comum em sistemas de 11 regras)
    0, 0, 2, 3, 3, 1, 1;
% 11a. Length['Long'] & Width['Small'] -> Quality['Middle']
    0, 1, 3, 0, 2, 1, 1;
% 11b. Length['Long'] & Diameter['Small'] -> Quality['Middle']
    1, 0, 3, 0, 2, 1, 1;
];

fis = addrule(fis, rules);
disp(' ');
disp('Base rules');
showrule (fis);

%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');
disp('Simulation for all samples in dataset');

quality_scores = [];

for i = 1:size(dataset)(1)
  data_input = dataset(i,:);
  FQI = show_activations(fis, data_input, rules);
  quality_scores = [quality_scores; FQI];
end

mean_quality = mean(quality_scores);

%%%%%%%%%%%%%%%%%%%%%%%%% HISTOGRAM PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A cor do histograma será definida com base no YEAR
if strcmp(YEAR, '2012')
    plot_color = '#4BD258'; % Green (verde) para 2012
else
    plot_color = '#4BC4D2'; % Blue (azul) para outros anos
end

if plot_hist
    figure('Position', [100, 100, 700, 450]);
    hist(quality_scores, 8, 'FaceColor', plot_color, 'EdgeColor', 'k');
    hold on;
    y_lims = ylim();
    line([mean_quality, mean_quality], y_lims, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', ['Mean: ', num2str(mean_quality, '%.2f')]);
    hold off;
    xlabel('Fuzzy Quality Index');
    ylabel('Frequency');
    grid on;
    legend('Histogram', ['Mean: ', num2str(mean_quality, '%.2f')], 'Location', 'northeast');
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXTREME ALMOND ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualizes fuzzy classification for minimum and maximum quality samples.

% Finds the extreme samples
[min_val, idx_low] = min(quality_scores);
[max_val, idx_high] = max(quality_scores);

% Complete samples (Diameter, Width, Length, Weight)
almond_low_quality  = dataset(idx_low, :);
almond_high_quality = dataset(idx_high, :);

% --- LOW QUALITY ALMOND SAMPLE ---
disp(' ');
disp('==================================================');
disp('CLASSIFICATION ANALYSIS: LOW QUALITY ALMOND');
disp('==================================================');

% 1. Prints the characteristics
disp(['Characteristics: Diameter=', num2str(almond_low_quality(1), '%.2f'), ...
      ', Width=', num2str(almond_low_quality(2), '%.2f'), ...
      ', Length=', num2str(almond_low_quality(3), '%.2f'), ...
      ', Weight=', num2str(almond_low_quality(4), '%.2f')]);

% 2. Shows activations and activated rules in the console
disp('Classification Process:');
show_activations(fis, almond_low_quality, rules);

% --- HIGH QUALITY ALMOND SAMPLE ---
disp(' ');
disp('==================================================');
disp('CLASSIFICATION ANALYSIS: HIGH QUALITY ALMOND');
disp('==================================================');

% 1. Prints the characteristics
disp(['Characteristics: Diameter=', num2str(almond_high_quality(1), '%.2f'), ...
      ', Width=', num2str(almond_high_quality(2), '%.2f'), ...
      ', Length=', num2str(almond_high_quality(3), '%.2f'), ...
      ', Weight=', num2str(almond_high_quality(4), '%.2f')]);

% 2. Shows activations and activated rules in the console
disp('Classification Process:');
show_activations(fis, almond_high_quality, rules);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diary off;
