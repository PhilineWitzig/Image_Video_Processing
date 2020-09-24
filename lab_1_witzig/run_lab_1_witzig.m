% Lab 1 - Philine Witzig and 23.09.2020


% getting path to project as it might vary from MATLAB user_path
path = matlab.desktop.editor.getActiveFilename;
path_split = strsplit(path, '/');
path_cur_folder = char(join(path_split(1:end - 1), '/'));

% choose the exercise to be executed
exercise = input('Enter the number of the exercise you want to execute: ');

switch exercise
    case 1
        %% Exercise 1 - "Images and Color Tables"

        % load images
        tree_file = fullfile(path_cur_folder, 'Images', 'trees.tif');
        lena_file = fullfile(path_cur_folder, 'Images', 'lena.tif');
        [I_trees, map_trees] = imread(tree_file);
        I_lena = imread(lena_file);

        % display original images
        disp('Original images');

        imshow(I_trees, map_trees);
        title('Trees original');
        pause();

        imshow(I_lena);
        title('Lena original');
        pause();

        % display images in gray scale
        disp('Images in grayscale');
        map_trees_gray = rgb2gray(map_trees);
        I_lena_gray = rgb2gray(I_lena);

        imshow(I_trees, map_trees_gray);
        title('Trees grayscale');
        pause();

        imshow(I_lena_gray);
        title('Lena grayscale');
        pause();

        % invert images
        disp('Inverted images')
        map_trees_gray_inv = imcomplement(map_trees_gray);
        I_lena_gray_inv = imcomplement(I_lena_gray);

        imshow(I_trees, map_trees_gray_inv);
        title('Trees grayscale inverted');
        pause();

        imshow(I_lena_gray_inv);
        title('Lena grayscale inverted');
        pause();

        % gamma correction
        disp('Gamma corrected images');
        gamma_values = [0.5 0.75 1.0 1.25 1.5 1.75 2.0]; % gamma test candidates

        for i=1:length(gamma_values)
            map_trees_gamma = gamma_corr(map_trees, gamma_values(i));
            imshow(I_trees, map_trees_gamma);
            title(strcat('Trees gamma corrected, gamma= ', num2str(gamma_values(i))));

            pause();
        end

        % create chess board
        I_chess = repmat([1 2; 2 1], 4, 4);
        cmap = [1 1 0; 0 0 1]; % [yellow, blue]
        map_chess = colormap(cmap);

        path_chess_indexed = fullfile(path_cur_folder, 'Images', 'chess_indexd.tif');
        path_chess_rgb = fullfile(path_cur_folder, 'Images', 'chess_rgb.tif');
        imwrite(I_chess, map_chess, path_chess_indexed);
        imwrite(ind2rgb(I_chess, map_chess), path_chess_rgb); % convert ind to rgb and save it

        
    case 2
        %% Exercise 2 - "Image Qunatization"

        lena_y_file = fullfile(path_cur_folder, 'Images', 'lena-y.png');
        I_lena_y = imread(lena_y_file);
        levels = [128 64 32 16 8 4 2];
        disp(length(unique(I_lena_y)));

        for i=1:length(levels)
            bin_size = 256 / (levels(i) - 1);
            I_lena_y_quant = bin_size * floor(I_lena_y ./ (bin_size)); % elementwise division
            cur_path = fullfile(path_cur_folder, 'Images', strcat('lena-y_', num2str(levels(i)), '.png'));
            imwrite(I_lena_y_quant, cur_path);
        end
        
        
    case 3
        %% Exercise 3 - "Filtering"
        
        I_gold_text = imread(fullfile(path_cur_folder, 'Images', 'gold-text.png'));
        
        F1 = [0.0357; 0.2411; 0.4464; 0.2411; 0.0357] * [0.0357 0.2411 0.4464 0.2411 0.0357];
        freqz2(F1);
        title('Frequency response of 5x5 Filter');
        pause();
        
        I_gold_text_filt = imfilter(I_gold_text, F1, 'conv');
        imshow(I_gold_text_filt);
        title('Test image convolved with 5x5 filter, edges trimmed');
        pause();
        
        F2 = (1 / 6) .* [-1 -4 -1; -4 26 -4; -1 -4 -1];
        freqz2(F2);
        title('Frequence response of 3x3 Filter');
        pause();
        I_gold_text_filt = imfilter(I_gold_text, F2, 'conv');
        imshow(I_gold_text_filt);
        title('Test image convolved with 3x3 filter, edges trimmed');
        pause();

        
    case 4
        %% Exercise 4 - "Correlation"
        I_gold_text = imread(fullfile(path_cur_folder, 'Images', 'gold-text.png'));
        height = size(I_gold_text, 1);
        width = size(I_gold_text, 2);
        template = imread(fullfile(path_cur_folder, 'Images', 'g-letter.png'));
        height_p = size(template, 1);
        width_p = size(template, 2);
        phi = zeros(size(I_gold_text));
        disp(height_p);
        disp(width_p);
        
        % TODO: apply zero padding to image
        for l = 1:height 
            for k = 1:width
                for l_p = floor(-(height_p / 2)):floor((height_p / 2))
                    for k_p = floor(-(width_p / 2)):floor((width_p / 2))
                        phi(l, k) = phi(l, k) + template(l_p, k_p) * I_gold_text(l + l_p, k + k_p);
                    end
                end
            end
        end
        
        imshow(phi, []);
        pause();

    case 5
        %% Exercise 5 - "Resampling"

    case 6
        %% Exercise 6 - "Phase and Magnitude of th 2DFT"

    case 7
        %% Exercise 7 - "Weber Law"

    otherwise
        disp('Invalid exercise number.')
        
end

%% Helper functions
function pause()
    try
        waitforbuttonpress;
    catch
        disp('Image closed.');
    end
end

function new_color_table = gamma_corr(color_table, gamma)
    new_color_table = color_table.^gamma; % elementwise exponential, gamma is the same for all channels
end
