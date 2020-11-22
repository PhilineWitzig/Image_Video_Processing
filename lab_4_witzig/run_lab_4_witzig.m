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
    T = cell2table(cell(0, 8));
    T.Properties.VariableNames = {'filename', 'point2point_MSE', 'point2point_HD', 'angular', 'rgb_1', 'rgb_2', 'yuv', 'br'};
    data = load("data/data.mat").data;

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
                    
                    % Get points of reference pc and pc under evaluation
                    ref = pc_ref.Location;
                    eval = pc_eval.Location;
                    
                    % Get nearest neighbors in both directions 
                    [idcs_re, dist_re] = knnsearch(ref, eval, 'Distance', 'euclidean'); % from  eval to ref
                    [idcs_er, dist_er] = knnsearch(eval, ref, 'Distance', 'euclidean'); % from  ref to eval
                    
                    % If normals don't exist yet, uncomment this code to
                    % compute normals (this may take a while)
                    % normals_ref = pcnormals(pc_ref, 128);
                    % pc_ref.Normal = normals_ref;
                    % normals_eval = pcnormals(pc_eval, 128);
                    % pc_eval.Normal = normals_eval;
                    % pcwrite(pc_ref, cur_ref_file_name);
                    % pcwrite(pc_eval, cur_eval_file_name);
                    
                    qual_point2point = Point2Point(dist_re, dist_er);
                    qual_plane2plane = Plane2Plane(pc_ref, pc_eval, idcs_re, idcs_er);
                    qual_color = ColorMetric(pc_ref, pc_eval, idcs_re, idcs_er);
                    bitrate = data(strcmp(data.filename, cur_eval_file_name), :).bpp;
                    T = [T;{cur_eval_file_name, qual_point2point(1), qual_point2point(2), qual_plane2plane, qual_color(1), qual_color(2), qual_color(1,3), bitrate}];
                    
                    writetable(T, "objective.csv");
                    
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

function qual = Point2Point(dist_re, dist_er)
    qual_MSE = max(sum(dist_re .^2) / length(dist_re), sum(dist_er .^2) / length(dist_er));
       
    hausdorff_er = max(min(dist_re,[],2)); % Directed from eval to ref
    hausdorff_re = max(min(dist_er, [], 2)); % Directed from ref to eval
    qual_HD = max(hausdorff_er, hausdorff_re);
    
    qual = [qual_MSE, qual_HD];
    
end

function qual = Plane2Plane(pc_ref, pc_eval, icds_re, idcs_er)
    
    normals_ref = pc_ref.Normal;
    normals_eval = pc_eval.Normal;

    % from eval to ref
    normals_ref_matched = normals_ref(icds_re, :);
    theta_hat = acos(dot(normals_ref_matched, normals_eval, 2) ./ (norm(normals_ref_matched) .* norm(normals_eval)));
    theta = min(theta_hat, pi - theta_hat);
    qual_re = 1 - (2 .* theta ./ pi);
    qual_re_MSE = sum(qual_re .^2) / length(qual_re);
    
    % from ref to eval
    normals_eval_matched = normals_eval(idcs_er, :);
    theta_hat = acos(dot(normals_eval_matched, normals_ref, 2) ./ (norm(normals_eval_matched) .* norm(normals_ref)));
    theta = min(theta_hat, pi - theta_hat);
    qual_er = 1 - (2 .* theta ./ pi);
    qual_er_MSE = sum(qual_er .^2) / length(qual_er);
    
    % get max
    qual = max(qual_re_MSE, qual_er_MSE);
end

function qual = ColorMetric(pc_ref, pc_eval, idcs_re, idcs_er)

    
    colors_ref_matched = double(pc_ref.Color(idcs_re, :));
    colors_eval_matched = double(pc_eval.Color(idcs_er, :));
    
    % PSNR RGB
    color_diff_ref = colors_ref_matched - double(pc_eval.Color);
    color_diff_eval = colors_eval_matched - double(pc_ref.Color);
    
    PSNR_RGB_1_ref = get_PSNR(color_diff_ref);
    PSNR_RGB_1_eval = get_PSNR(color_diff_eval);
    PSNR_RGB_1 = max(PSNR_RGB_1_ref, PSNR_RGB_1_eval);
    
    PSNR_RGB_2_ref = (get_PSNR(color_diff_ref(:, 1)) + get_PSNR(color_diff_ref(:, 2)) + get_PSNR(color_diff_ref(:, 3))) / 3;
    PSNR_RGB_2_eval = (get_PSNR(color_diff_eval(:, 1)) + get_PSNR(color_diff_eval(:, 2)) + get_PSNR(color_diff_eval(:, 3))) / 3;
    PSNR_RGB_2 = max(PSNR_RGB_2_ref, PSNR_RGB_2_eval);
    
    % PSNR YUV
    color_diff_ref_YUV = RGB2YUV(colors_ref_matched) - RGB2YUV(double(pc_eval.Color));
    PSNR_YUV_ref = (6 * get_PSNR(color_diff_ref_YUV(:, 1)) + get_PSNR(color_diff_ref_YUV(:, 2)) + get_PSNR(color_diff_ref_YUV(:, 3))) / 8;
    color_diff_eval_YUV = RGB2YUV(colors_eval_matched) - RGB2YUV(double(pc_ref.Color));
    PSNR_YUV_eval = (6 * get_PSNR(color_diff_eval_YUV(:, 1)) + get_PSNR(color_diff_eval_YUV(:, 2)) + get_PSNR(color_diff_eval_YUV(:, 3))) / 8;
    PSNR_YUV = max(PSNR_YUV_ref, PSNR_YUV_eval);
    
    qual = [PSNR_RGB_1, PSNR_RGB_2, PSNR_YUV];
    
end

function colors_YUV = RGB2YUV(colors_RGB)
    colors_YUV = zeros(size(colors_RGB));
    colors_YUV(:, 1) = 0.299 .* colors_RGB(:, 1) + 0.587 .* colors_RGB(:, 2) + 0.144 .* colors_RGB(:, 3);
    colors_YUV(:, 2) = 0.493 .* (colors_RGB(:, 3) - colors_YUV(:, 1));
    colors_YUV(:, 3) = 0.877 .* (colors_RGB(:, 1) - colors_YUV(:, 1));
end

function PSNR = get_PSNR(colors)

    [height, width] = size(colors);
    MSE = sum(colors.^2, 'all') / (height * width);
    PSNR = 10 * log10(255^2 / MSE);
end