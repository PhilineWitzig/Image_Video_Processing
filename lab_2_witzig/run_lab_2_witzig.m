% Lab 2 - Philine Witzig 14.10.2020

% getting path to project as it might vary from MATLAB user_path
path = matlab.desktop.editor.getActiveFilename;
path_split = strsplit(path, '/');
path_cur_folder = char(join(path_split(1:end - 1), '/'));


% choose the exercise to be executed
exercise = input("Enter the number of the exercise you want to execute: ");

switch exercise
    case 1
    %% Exercise 1 - "Fixed Threshold Method"
    % Lena image
    I = imread("Images/lena-y.png");
    % compute the dynamic range of the image, half it and normalize to
    % range [0, 1]
    t = double((max(I, [], 'all') - min(I, [], 'all')) / 2) / 255;
    BW = im2bw(I, t);
    figure('name', strcat("Fixed Thresh ", num2str(t)));
    imshow(BW);
    
    % Wool image
    I = imread("Images/wool.png");
    t = double((max(I, [], 'all') - min(I, [], 'all')) / 2) / 255;
    BW = im2bw(I, t);
    figure('name', strcat("Fixed Thresh ", num2str(t)));
    imshow(BW);

    case 2
    %% Exercise 2 - "Random Threshold Method"
    % Lena image
    I = imread("Images/lena-y.png");
    [height, width] = size(I);
    noise = 0.2 * unidrnd(255, height, width);
    I_n = mat2gray(double(I)+ noise);
    t = (max(I_n, [], 'all') - min(I_n, [], 'all')) / 2;
    BW = im2bw(I_n, t);
    figure('name', "Random Thresh Method");
    imshow(BW);
    
    % Wool image
    I = imread("Images/wool.png");
    [height, width] = size(I)
    noise = 0.6 * unidrnd(255, height, width);
    I_n = mat2gray(double(I)+ noise);
    t = (max(I_n, [], 'all') - min(I_n, [], 'all')) / 2;
    BW = im2bw(I_n, t);
    figure('name', "Random Thresh Method");
    imshow(BW);

    case 3
    %% Exercise 3 "Ordered Threshold Method"
    I = imread("Images/lena-y.png");
    I_quant = quantize(double(I), 0, 5);
    M = [1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5];
    I_thresh = ordered_thresh(I_quant, M);
    figure('name', "Test filtering");
    imshow(I_thresh);

    case 4
    %% Exercise 4 - "Ordered Matrix with Centered Points"
    % lena image
    I = imread("Images/lena-y.png");
    C6 = [34 25 21 17 29 33; 30 13 9 5 12 24; 18 6 1 0 8 20; 22 10 2 3 4 16; 26 14 7 11 15 28; 35 31 19 23 27 32];
    E6 = [30 22 16 21 33 35; 24 11 7 9 26 28; 13 5 0 2 14 19; 15 3 1 4 12 18; 27 8 6 10 25 29; 32 20 17 23 31 34];

    N_c6 = length(unique(C6));
    N_e6 = length(unique(E6));
    
    I_quant_c6 = quantize(double(I), 0, N_c6);
    I_quant_e6 = quantize(double(I), 0, N_e6);
    
    I_thresh_c6 = ordered_thresh(I_quant_c6, C6);
    I_thresh_e6 = ordered_thresh(I_quant_e6, E6);
    
    figure('name', "Ordered matrix centered points C6");
    imshow(I_thresh_c6);
    
    figure('name', "Ordered matrix centered points E6");
    imshow(I_thresh_e6);
   
    % wool img 
    I = imread("Images/wool.png");
    
    I_quant_c6 = quantize(double(I), 0, N_c6);
    I_quant_e6 = quantize(double(I), 0, N_e6);
    
    I_thresh_c6 = ordered_thresh(I_quant_c6, C6);
    I_thresh_e6 = ordered_thresh(I_quant_e6, E6);
    
    figure('name', "Ordered matrix centered points C6");
    imshow(I_thresh_c6);
    
    figure('name', "Ordered matrix centered points E6");
    imshow(I_thresh_e6);
    
    case 5
    %% Exercise 5 - "Diagonal Ordered Matrix with Balanced Centered Points"
    O81 = [13 9 5 12; 6 1 0 8; 10 2 3 4; 14 7 11 15];
    O82 = [18 22 26 19; 25 30 31 23; 21 29 28 27; 17 24 20 16];
    O8 = [O81 O82; O82 O81];
    N_O8 = length(unique(O8));
    
    I = imread("Images/lena-y.png");
    I_quant_O8 = quantize(I, 0, N_O8);
    I_thresh_O8 = ordered_thresh(I_quant_O8, O8);
    figure('name', "Diagonal ordered matrix with balanced cenetered points O8");
    imshow(I_thresh_O8);
    
    I = imread("Images/wool.png");
    I_quant_O8 = quantize(I, 0, N_O8);
    I_thresh_O8 = ordered_thresh(I_quant_O8, O8);
    figure('name', "Diagonal ordered matrix with balanced cenetered points O8");
    imshow(I_thresh_O8);
    
    
    case 6
    %% Exercise 6 - "Ordered Matrix with Dispersed Dots"

    case 7
    %% Exercise 7 - "Error diffusion Method"

end

function I_thresh = ordered_thresh(I, M)
    [height, width] = size(I);
    [height_m, width_m] = size(M);
    I_thresh = zeros(height, width);
    for x=0:width_m:(width - width_m - 1)
        for y=0:height_m:(height - height_m - 1)
            for x_m=1:width_m
                for y_m=1:height_m
                    I_thresh(y + y_m, x + x_m) =  (I(y + y_m, x + x_m) > M(y_m, x_m)) * 1 + (~(I(y + y_m, x + x_m) > M(y_m, x_m))) * 0;
                end
            end
        end
    end
end

function I_quant = quantize(I, lower, upper)
    I_min = min(I, [], 'all');
    I_max = max(I, [], 'all');
    I_quant = uint8((((I - I_min)) .* (upper - lower)) ./ ((I_max - I_min) + lower));
end