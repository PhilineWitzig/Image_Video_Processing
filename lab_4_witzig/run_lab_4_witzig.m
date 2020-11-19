% Lab 4 - Philine Witzig 17.11.2020

%pc = pcread("data/");
%pcshow(pc);


subj_data = readmatrix("data/subjective_results.csv");
[height, width] = size(subj_data);
% drop model names and sciper numbers
subj_data = subj_data(2:height, 2:width);

exercise = input("Enter the number of the exercise you want to execute: ");

switch exercise
    case 1
    %% Exercise 2.3 - "Subjective Quality Assessment"
    subj_data_cleaned = detect_outliers(subj_data);
    while (subj_data_cleaned ~= subj_data)
        subj_data_cleaned = detect_outliers(subj_data);
    end 

    % Compute DMOS and CI
    [DMOS, DV] = get_DMOS(subj_data_cleaned);
    CI = get_CI(DV);
    % Create plots
    plot_DMOS(subj_data_cleaned, DMOS, CI);

    
    case 2
    %% Exercise 2.4 - "Objective Quality Assessment"
    myDir = pwd;
    pc_files = dir(fullfile(myDir, 'models')); % get all files

    for f_ref=1:length(pc_files)
        cur_ref_file_name = pc_files(f_ref).name;
        
        if contains(cur_ref_file_name, "ref")
            pc_ref = pcread(fullfile("models", cur_ref_file_name));
            cur_model = strsplit(cur_ref_file_name, "_");
            cur_model = cur_model(1);
            
            for f_eval=1:length(pc_files)
                cur_eval_file_name = pc_files(f_eval).name;
                if contains(cur_eval_file_name, cur_model) && ~(contains(cur_eval_file_name, "ref"))
                    pc_eval = pcread(fullfile("models", cur_eval_file_name));
                    disp("Distance between: ")
                    disp(cur_ref_file_name);
                    disp(cur_eval_file_name);
                    quality = symm_P2P(pc_ref, pc_eval, ["MSE", "Hausdorff"]);
                    disp(strcat("Objective quality for MSE, Huassdorff: ", num2str(quality))); 
                    disp("_____");
                end
            end
        end
    end
        
    case 3
    %% Exercise 2.5 - "Benchmarking of Objective Quality Metrics"

end
    
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
    DMOS = zeros(stimuli, 1); 
    DV = zeros(stimuli, participants);
    for r=1:9:stimuli
        V_REF = T(r, :);
        for d=1:8
            DV(r + d, :) = (T(r + d, :) - V_REF) + 5;
            DMOS(r + d) = sum(DV(r + d, :)) / participants;
        end
    end

    % without the reference stimuli
    for i=1:8:32
        DMOS(i, :) = [];
        DV(i, :) = [];
    end
end

function CI=get_CI(DV)
    [stimuli, participants] = size(DV);
    CI = zeros(stimuli, 1); % value of reference value is not included
    for j=1:stimuli
        CI(j, :) = icdf('T', 1 - (0.95/2), participants - 1) .* (std(DV(j, :)) / sqrt(participants));
    end
end

function plot_DMOS(T, DMOS, CI)
   
    [stimuli, participants] = size(DMOS);
    
    % remove reference from DMOS and CI
    
    data = load("data/data.mat");
    whos data
    disp(data.data)
    
    for i=1:8:stimuli
        codec1 = DMOS(i:i+3, :);
        codec2 = DMOS(i+4:i+7, :);
        x_codec1 = data.data.bpp(i:i+3);
        x_codec2 = data.data.bpp(i+4:i+7);
        if i == 1
            figure('name', "Long dress");
        elseif i == 9
            figure('name', "Guanyin");
        elseif i == 17
            figure('name', "Phil");   
        else
            figure('name', "Rhetorician");
        end
        errorbar(x_codec1, codec1, CI(i:i+3));
        hold on;
        errorbar(x_codec2, codec2, CI(i+4:i+7), 'r');
        ylim([0 6]);
        yticks([1 2 3 4 5]);
        xlabel("Bitrates");
        ylabel("DMOS");
        legend({"Codec 1: pcc geo color", "Codec 2: cwi pcl"}, 'location', 'northwest');
    end
end

function qual=symm_P2P(pc1, pc2, metrics)
    ref = pc1.Location;
    eval = pc2.Location;
    qual=[];
    [idcs, dist] = knnsearch(ref, eval, 'Distance', 'euclidean');
  
    for metric = metrics
        if metric == "MSE"
            qual=[qual, sum(dist .^2) / length(dist)];
        elseif metric == "Hausdorff"
            qual=[qual, max(min(dist,[],2))]; % Directed from eval to ref
            % hba = max(min(D));% Directed from ref to eval
            % H = max([hab,hba]);
        else
            disp("Invalid metric")
            return
        end
    end
    
end