% Lab 1 - Philine Witzig and 23.09.2020

%% Exercise 1 - "Images and Color Tables"

% getting path to project as it might vary from MATLAB user_path
path = matlab.desktop.editor.getActiveFilename;
path_split = strsplit(path, '/');
path_cur_folder = char(join(path_split(1:end - 1), '/'));

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

% Helper functions
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


%% Exercise 2 - "Image Qunatization"



%% Exercise 3 - "Filtering"


%% Exercise 4 - "Correlation"


%% Exercise 5 - "Resampling"


%% Exercise 6 - "Phase and Magnitude of th 2DFT"


%% Exercise 7 - "Weber Law"