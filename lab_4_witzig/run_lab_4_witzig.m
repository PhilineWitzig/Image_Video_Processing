% Lab 4 - Philine Witzig 17.11.2020
addpath('data/');
import angularSimilarity.*

subj_data = readtable("data/subjective_results.csv");
subj_data(1, :) = [];

exercise = input("Enter the number of the exercise you want to execute: ");
subj_data_cleaned = detect_outliers(subj_data);

while (~isempty(setdiff(subj_data_cleaned,subj_data)))
    subj_data = subj_data_cleaned;
    subj_data_cleaned = detect_outliers(subj_data);
end 

switch exercise
    case 1
    %% Exercise 2.3 - "Subjective Quality Assessment"

    % Compute DMOS and CI
    [DMOS, DV] = get_DMOS(subj_data_cleaned);
    CI = get_CI(DV);
    % Create plots
    plot_DMOS(DMOS, CI);

    case 2
    %% Exercise 2.4 - "Objective Quality Assessment"

    if isfile("objective.csv")
        T = readtable("objective.csv");
    else
        myDir = pwd;
        pc_files = dir(fullfile(myDir, 'models')); % get all files
        % store the objective metric values in a table
        T = cell2table(cell(0, 8));
        T.Properties.VariableNames = {'filename', 'point2point_MSE', 'point2point_HD', 'angular', 'rgb_1', 'rgb_2', 'yuv', 'br'};
        data = load("data/data.mat").data;
        for f_ref=1:length(pc_files)
            % get reference file
            cur_ref_file_name = pc_files(f_ref).name;

            if contains(cur_ref_file_name, "ref")
                pc_ref = pcread(fullfile("models", cur_ref_file_name));
                cur_model = strsplit(cur_ref_file_name, "_");
                cur_model = cur_model(1);

                for f_eval=1:length(pc_files)
                    % get evaluation file
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
                        % compute normals (this may take a while, thus we have them precomputed)
                        % normals_ref = pcnormals(pc_ref, 128);
                        % pc_ref.Normal = normals_ref;
                        % normals_eval = pcnormals(pc_eval, 128);
                        % pc_eval.Normal = normals_eval;
                        % pcwrite(pc_ref, cur_ref_file_name);
                        % pcwrite(pc_eval, cur_eval_file_name);

                        % apply the different metrics
                        qual_point2point = Point2Point(dist_re, dist_er);
                        qual_plane2plane = Plane2Plane(pc_ref, pc_eval, idcs_re, idcs_er);
                        qual_color = ColorMetric(pc_ref, pc_eval, idcs_re, idcs_er);
                        bitrate = data(strcmp(data.filename, cur_eval_file_name), :).bpp;
                        
                        % store results in table
                        T = [T;{cur_eval_file_name, qual_point2point(1), qual_point2point(2), qual_plane2plane, qual_color(1), qual_color(2), qual_color(1,3), bitrate}];
  
                    end
                end
            end
        end
        % store table
        writetable(T, "objective.csv");
    end
    disp(T);
    plot_objectives(T);

    case 3
    %% Exercise 2.5 - "Benchmarking of Objective Quality Metrics"
    plot_DMOS_objective("dmos.csv", "objective.csv");

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
    MOS_r2_y = zeros(n_r2, participants-1);
 
    for model=1:stimuli
        MOS_r1(model) = mean(table2array(T(model, 2:participants)));
    end
    
    for content=1:9:stimuli
        MOS_r2_x(content) = mean(table2array(T(content:content + 8, 2:participants)), 'all');
        MOS_r2_y(content, :) = mean(table2array(T(content:content + 8, 2:participants)));
    end
    
    % Init PLCC array
    PLCC_r1 = zeros(participants, 1);
    PLCC_r2 = zeros(participants, 1);
    
    for subj=2:participants
        score_subj = table2array(T(:, subj));
        nominator_r1 = n_r1 * sum(MOS_r1.* score_subj) - sum(MOS_r1) * sum(score_subj);
        denominator_r1 = sqrt(n_r1 * sum(MOS_r1.^2) - sum(MOS_r1)^2) * sqrt(n_r1 * sum(score_subj.^2) - sum(score_subj)^2);
        PLCC_r1(subj) = nominator_r1/denominator_r1;
        
        % compute MOS
        MOS_subj = MOS_r2_y(:, subj-1);
        nominator_r2 = n_r2 * sum(MOS_r2_x .* MOS_subj) - sum(MOS_r2_x) * sum(MOS_subj);
        denominator_r2 = sqrt(n_r2 * sum(MOS_r2_x.^2) - sum(MOS_r2_x)^2) * sqrt(n_r2 * sum(MOS_subj.^2) - sum(MOS_subj)^2);
        PLCC_r2(subj) = nominator_r2/denominator_r2;
    end
    
    % Remove the worst outlier
    worst_value = 0.0;
    worst_candidate = 0;
    for i=2:participants
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
    % Compute DMOS and DV
    [stimuli, participants] = size(T);
    DMOS = zeros(stimuli, 1); 
    DV = zeros(stimuli, participants - 1);
    for r=1:9:stimuli
        V_REF = table2array(T(r, 2:participants));
        for d=1:8
            DV(r + d, :) = (table2array(T(r + d, 2:participants)) - V_REF) + 5;
            DMOS(r + d) = sum(DV(r + d, :)) / (participants - 1);
        end
    end

    % without the reference stimuli
    for i=1:8:32
        DMOS(i, :) = [];
        DV(i, :) = [];
    end
end

function CI=get_CI(DV)
    % Compute the confidence intervals, alpha=0.05
    [stimuli, participants] = size(DV);
    CI = zeros(stimuli, 1); % value of reference value is not included
    for j=1:stimuli
        CI(j, :) = icdf('T', 1 - (.05/2), participants - 1) .* (std(DV(j, :)) / sqrt(participants));
    end
end

function plot_DMOS(DMOS, CI)
    % Produce the DMOS vs bitrate plots
    
    [stimuli, participants] = size(DMOS);
    
    % remove reference from DMOS and CI
    data = load("data/data.mat");
    whos data
    disp(data.data)
    
    % store data
    names = data.data.filename;
    T = array2table([names,  DMOS, CI, data.data.bpp], 'VariableNames',  {'filename','dmos', 'ci', 'br'});
    writetable(T, "dmos.csv");
    
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
        % swap ordering here to make color coding consistent with other
        % tasks
        errorbar(x_codec2, codec2, CI(i+4:i+7));
        hold on;
        errorbar(x_codec1, codec1, CI(i:i+3), 'r');
        ylim([0 6]);
        yticks([1 2 3 4 5]);
        xlabel("Bitrates");
        ylabel("DMOS");
        legend({"Codec 1: cwi pcl", "Codec 2: pcc geo color"}, 'location', 'northwest');
    end
end

function qual = Point2Point(dist_re, dist_er)
    % Get the point 2 point metric values using MSE and Hausdorff
    qual_MSE = max(sum(dist_re .^2) / length(dist_re), sum(dist_er .^2) / length(dist_er));
       
    hausdorff_er = max(min(dist_re,[],2)); % Directed from eval to ref
    hausdorff_re = max(min(dist_er, [], 2)); % Directed from ref to eval
    qual_HD = max(hausdorff_er, hausdorff_re);
    
    qual = [qual_MSE, qual_HD];
    
end

function qual = Plane2Plane(pc_ref, pc_eval)
    % Get the plane to plane metric value using angular similarity
    qual = angularSimilarity(pc_ref, pc_eval, 'MSE');

end

function qual = ColorMetric(pc_ref, pc_eval, idcs_re, idcs_er)
    % Computer the color metric values on RGB and YUV space
    % Note that there are two version for computing the PSNR on RGB
    
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
    % Get the YUV values from RGB
    colors_YUV = zeros(size(colors_RGB));
    colors_YUV(:, 1) = 0.299 .* colors_RGB(:, 1) + 0.587 .* colors_RGB(:, 2) + 0.144 .* colors_RGB(:, 3);
    colors_YUV(:, 2) = 0.493 .* (colors_RGB(:, 3) - colors_YUV(:, 1));
    colors_YUV(:, 3) = 0.877 .* (colors_RGB(:, 1) - colors_YUV(:, 1));
end

function PSNR = get_PSNR(colors)
    % Compute the PSNR according to the formula in the assignment
    [height, width] = size(colors);
    MSE = sum(colors.^2, 'all') / (height * width);
    PSNR = 10 * log10(255^2 / MSE);
end

function plot_objectives(T)
    % Produce plots for objective metric values
    
    data = table2array(T(:, 2:8));
    [row, column] = size(T);
    for metric=1:6
        for i=1:8:row
            if i == 1
                figure('name', "Guanyin");
            elseif i == 9
                figure('name', "Longdress");
            elseif i == 17
                figure('name', "Phil");   
            else
                figure('name', "Rhetorician");
            end
            codec1 = data(i:i+3,metric);
            br_1 = data(i:i+3, 7);
            codec2 = data(i+4:i+7, metric);
            br_2 = data(i+4:i+7, 7);
            plot(br_1, codec1);
            title(T.Properties.VariableNames{metric + 1});
            xlabel("Bitrates");
            ylabel("Objective Metric");
            hold on;
            plot(br_2, codec2, 'r');
            legend({"Codec 1: cwi pcl", "Codec 2: pcc geo color"}, 'location', 'northeastoutside');
        end
    end
end

function plot_DMOS_objective(dmos_file, obj_file)
    % Produce plots for DMOS vs objective metric values
    % Fit the objective values
    
    T_DMOS = readtable(dmos_file);
    T_OBJ = readtable(obj_file);
    T = join(T_OBJ, T_DMOS, 'Keys','filename');
    [models, metrics] = size(T_OBJ);
    
    T_corr = cell2table(cell(0, 10));
    T_corr.Properties.VariableNames = {'metric', 'pearson', 'spearman', 'rmse', 'pearson_fit_linear', 'spearman_fit_linear', 'rmse_fit_linear', 'pearson_fit_cubic', 'spearman_fit_cubic', 'rmse_fit_cubic'};
   
    for metric=2:metrics-1
        obj_scores = table2array(T(:, metric)); % metric values
        subj_scores = T.dmos;
        
        % produce plot
        figure;
        e = errorbar(obj_scores, subj_scores, T.ci);
        e.LineStyle = 'none';
        title(T.Properties.VariableNames{metric});
        xlabel("Objective metric");
        ylabel("DMOS");
        hold on;
      
        % get correlation coeefficients and RMSEs
        pearson = corr(subj_scores, obj_scores, 'Type', 'Pearson');
        spearman = corr(subj_scores, obj_scores, 'Type', 'Spearman');
        rmse = sqrt(sum((obj_scores - subj_scores).^2, 'all'));
        
        % fit data 
        linear_fit = polyfit(obj_scores, subj_scores, 1);
        cubic_fit = polyfit(obj_scores, subj_scores, 3);
        
        obj_scores_linear = polyval(linear_fit, obj_scores);
        obj_scores_cubic = polyval(cubic_fit, obj_scores);
        
        % reproduce plots on fitted data
        figure;
        e = errorbar(obj_scores_linear, T.dmos, T.ci);
        e.LineStyle = 'none';
        title(T.Properties.VariableNames{metric});
        xlabel("Objective metric, linear fit");
        ylabel("DMOS");
        hold on;
        
        figure;
        e = errorbar(obj_scores_cubic, T.dmos, T.ci);
        e.LineStyle = 'none';
        title(T.Properties.VariableNames{metric});
        xlabel("Objective metric, cubic fit");
        ylabel("DMOS");
        hold on;
        
        % recompute coefficients on fitted data
        pearson_linear = corr(subj_scores, obj_scores_linear, 'Type', 'Pearson');
        spearman_linear = corr(subj_scores, obj_scores_linear, 'Type', 'Spearman');
        rmse_linear = sqrt(sum((obj_scores_linear - subj_scores).^2, 'all'));
        
        pearson_cubic = corr(subj_scores, obj_scores_cubic, 'Type', 'Pearson');
        spearman_cubic = corr(subj_scores, obj_scores_cubic, 'Type', 'Spearman');
        rmse_cubic = sqrt(sum((obj_scores_cubic - subj_scores).^2, 'all'));
        
        % store results in table
        T_corr = [T_corr; {T.Properties.VariableNames{metric}, pearson, spearman, rmse, pearson_linear, spearman_linear, rmse_linear, pearson_cubic, spearman_cubic, rmse_cubic}];
        
    end
    
    writetable(T_corr, "correlation.csv");

end
    