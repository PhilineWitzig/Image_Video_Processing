% Lab 4 - Philine Witzig 17.11.2020

%pc = pcread("data/");
%pcshow(pc);

subj_data = readmatrix("data/subjective_results.csv");
[height, width] = size(subj_data);
% drop model names and sciper numbers
subj_data = subj_data(2:height, 2:width);

%% Exercise 2.3 - "Subjective Quality Assessment"

subj_data_cleaned = detect_outliers(subj_data);
while (subj_data_cleaned ~= subj_data)
    subj_data_cleaned = detect_outliers(subj_data);
end 

% Compute DMOS and CI
[DMOS, DV] = get_DMOS(subj_data_cleaned);
CI = get_CI(subj_data_cleaned, DV);
plot_DMOS(subj_data_cleaned, DMOS, CI);


%% Exercise 2.4 - "Objective Quality Assessment"


%% Exercise 2.5 - "Benchmarking of Objective Quality Metrics"


function T=detect_outliers(T)
    [stimuli, participants] = size(T);
    n_r1 = stimuli;
    n_r2 = 4;
    % MOS of all the subjects for a given stimulus i
    MOS_r1 = zeros(n_r1, 1);
    % MOS of all subjects for all the stimuli for a given content i
    MOS_r2_x = zeros(n_r2, 1);
    % MOS of one subject for all the stimuli for given content i
    MOS_r2_y = zeros(n_r2, participants);
 
    for model=1:stimuli
        MOS_r1(model) = mean(T(model, :));
    end
    
    for content=1:9:stimuli
        MOS_r2_x(content) = mean(T(content:content + 8, :), 'all');
        MOS_r2_y(content, :) = mean(T(content:content + 8, :));
    end
    
    % get r1
    PLCC_r1 = zeros(participants, 1);
    PLCC_r2 = zeros(participants, 1);
    
    for subj=1:participants
        score_subj = T(:, subj);
        nominator_r1 = n_r1 * sum(MOS_r1.* score_subj) - sum(MOS_r1) * sum(score_subj);
        denominator_r1 = sqrt(n_r1 * sum(MOS_r1.^2) - sum(MOS_r1)^2) * sqrt(n_r1 * sum(score_subj.^2) - sum(score_subj)^2);
        PLCC_r1(subj) = nominator_r1/denominator_r1;
        
        MOS_subj = MOS_r2_y(:, subj);
        nominator_r2 = n_r2 * sum(MOS_r2_x .* MOS_subj) - sum(MOS_r2_x) * sum(MOS_subj);
        denominator_r2 = sqrt(n_r2 * sum(MOS_r2_x.^2) - sum(MOS_r2_x)^2) * sqrt(n_r2 * sum(MOS_subj.^2) - sum(MOS_subj)^2);
        PLCC_r2(subj) = nominator_r2/denominator_r2;
    end
    
    worst_value = 0.0;
    worst_candidate = 0;
    for i=1:participants
        if PLCC_r1(i) < 0.75 && PLCC_r2(i) < 0.8
            diff = mean((0.75- PLCC_r1(i)) + (0.8 - PLCC_r2(i)));
            if diff > worst_value
                worst_value = diff;
                worst_candidate = i;
            end
        end
    end
    
   
    if worst_candidate ~= 0
        T(worst_candidate) = [];
    end
end


function [DMOS,DV]=get_DMOS(T)
    [stimuli, participants] = size(T);
    DMOS = zeros(stimuli, 1); % without the reference stimuli? Should be 0 here
    DV = zeros(stimuli, participants);
    for r=1:9:stimuli
        V_REF = T(r, :);
        for d=0:8
            DV(r + d, :) = (T(r + d, :) - V_REF) + 5;
            DMOS(r + d) = sum(DV(r + d, :)) / participants;
        end
    end
end

function CI=get_CI(T, DV)
    [stimuli, participants] = size(T);
    CI = zeros(stimuli, 1); % value of reference value is not 
    for j=1:stimuli
        CI(j, :) = icdf('T', 1 - (0.95/2), participants - 1) .* (std(DV(j, :)) / sqrt(participants));
    end
end


function plot_DMOS(T, DMOS, CI)
    [stimuli, participants] = size(T);
    for i=1:9:stimuli
        codec1 = DMOS(i+1:i+4, :);
        codec2 = DMOS(i+5:i+8, :);
        x = 1:4;
        figure('name', "Errorbars");
        errorbar(codec1, x, CI(i+1:i+4));
        hold on;
        errorbar(codec2, x, CI(i+5:i+8), 'r');
    end
end
