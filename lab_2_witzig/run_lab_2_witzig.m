% Lab 2 - Philine Witzig 14.10.2020

% getting path to project as it might vary from MATLAB user_path
path = matlab.desktop.editor.getActiveFilename;
path_split = strsplit(path, '/');
path_cur_folder = char(join(path_split(1:end - 1), '/'));


% choose the exercise to be executed
exercise = input("Enter the number of the exercise you want to execute: ");
I_lena = imread("Images/lena-y.png");
I_wool = imread("Images/wool.png");

switch exercise
    case 1
    %% Exercise 1 - "Fixed Threshold Method"
    % Lena image
    % compute the dynamic range of the image, half it and normalize to
    % range [0, 1]
    t = double((max(I_lena, [], 'all') - min(I_lena, [], 'all')) / 2) / 255;
    BW = imbinarize(I_lena, t);
    figure('name', strcat("Fixed Thresh ", num2str(t)));
    imshow(BW);
    
    % Wool image
    t = double((max(I_wool, [], 'all') - min(I_wool, [], 'all')) / 2) / 255;
    BW = imbinarize(I_wool, t);
    figure('name', strcat("Fixed Thresh ", num2str(t)));
    imshow(BW);

    case 2
    %% Exercise 2 - "Random Threshold Method"
    % Lena image
    [height, width] = size(I_lena);
    noise = 0.2 * unidrnd(255, height, width);
    I_n = mat2gray(double(I_lena)+ noise);
    t = (max(I_n, [], 'all') - min(I_n, [], 'all')) / 2;
    BW = im2bw(I_n, t);
    figure('name', "Random Thresh Method");
    imshow(BW);
    
    % Wool image
    [height, width] = size(I_wool);
    noise = 0.6 * unidrnd(255, height, width);
    I_n = mat2gray(double(I_wool)+ noise);
    t = (max(I_n, [], 'all') - min(I_n, [], 'all')) / 2;
    BW = imbinarize(I_n, t);
    figure('name', "Random Thresh Method");
    imshow(BW);

    case 3
    %% Exercise 3 "Ordered Threshold Method"
    % TODO remove this part actually, since we only had to implement the
    % functions here
    I_quant = quantize(double(I_lena), 0, 5);
    M = [1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5];
    I_thresh = ordered_thresh(I_quant, M);
    figure('name', "Test filtering");
    imshow(I_thresh);

    case 4
    %% Exercise 4 - "Ordered Matrix with Centered Points"
    C6 = [34 25 21 17 29 33; 30 13 9 5 12 24; 18 6 1 0 8 20; 22 10 2 3 4 16; 26 14 7 11 15 28; 35 31 19 23 27 32];
    E6 = [30 22 16 21 33 35; 24 11 7 9 26 28; 13 5 0 2 14 19; 15 3 1 4 12 18; 27 8 6 10 25 29; 32 20 17 23 31 34];

    N_c6 = length(unique(C6));
    N_e6 = length(unique(E6));
    
    % lena img
    I_quant_c6 = quantize(double(I_lena), 0, N_c6 - 1);
    I_quant_e6 = quantize(double(I_lena), 0, N_e6 - 1);
    
    I_thresh_c6 = ordered_thresh(I_quant_c6, C6);
    I_thresh_e6 = ordered_thresh(I_quant_e6, E6);
    
    figure('name', "Ordered matrix centered points C6");
    imshow(I_thresh_c6);
    
    figure('name', "Ordered matrix centered points E6");
    imshow(I_thresh_e6);
   
    % wool img 
    I_quant_c6 = quantize(double(I_wool), 0, N_c6 - 1);
    I_quant_e6 = quantize(double(I_wool), 0, N_e6 - 1);
    
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
    
    % lena image
    I_quant_O8 = quantize(double(I_lena), 0, N_O8 - 1);
    I_thresh_O8 = ordered_thresh(I_quant_O8, O8);
    figure('name', "Diagonal ordered matrix with balanced cenetered points O8");
    imshow(I_thresh_O8);
    
    % wool image
    I_quant_O8 = quantize(double(I_wool), 0, N_O8 - 1);
    I_thresh_O8 = ordered_thresh(I_quant_O8, O8);
    figure('name', "Diagonal ordered matrix with balanced cenetered points O8");
    imshow(I_thresh_O8);
    
    
    case 6
    %% Exercise 6 - "Ordered Matrix with Dispersed Dots"

    D3 = [8 4 5; 3 0 1; 7 2 6];
    U3 = ones(3, 3);
    D6 = [4.*D3 4.*D3+2.*U3; 4.*D3+3.*U3 4.*D3+U3];
 
    N_D6 = length(unique(D6));
    
    % Lena image
    I_quant_D6 = quantize(double(I_lena), 0, N_D6 - 1);
    I_thresh_D6 = ordered_thresh(I_quant_D6, D6);
    figure('name', "Diagonal ordered matrix with balanced cenetered points D6");
    imshow(I_thresh_D6);
    
    % Wool image
    I_quant_D6 = quantize(double(I_wool), 0, N_D6 - 1);
    I_thresh_D6 = ordered_thresh(I_quant_D6, D6);
    figure('name', "Diagonal ordered matrix with balanced cenetered points D6");
    imshow(I_thresh_D6);
    
    case 7
    %% Exercise 7 - "Error diffusion Method"
    % Binary thresholding based on dynamic range of image
    t = double((max(I_lena, [], 'all') - min(I_lena, [], 'all')) / 2) / 255;
 
    % filter entries are of format y, x, coeff
    filter_floyd = [0 1 (7/16); 1 -1 (3/16); 1 0 (5/16); 1 1 (1/16)];
    filter_stucki = [0 1 (8/42); 0 2 (4/42); 1 -2 (2/42); 1 -1 (4/42); 1 0 (8/42); 1 1 (4/42); 1 2 (2/42); 2 -2 (1/42); 2 -1 (2/42); 2 0 (4/42); 2 1 (2/42); 2 2 (1/42)];
    
    % Floyd Filter
    I = rescale(I_lena);  % Rescale input image to range [0, 1]
    I = error_diffusion(I, filter_floyd, t);
    figure('name', "Error diffusion method with Floyd & Steinberg");
    imshow(I);
    
    % Stucki Filter
    I = rescale(I_lena);  % Rescale input image to range [0, 1]
    I = error_diffusion(I, filter_stucki, t);
    figure('name', "Error diffusion method with Stucki");
    imshow(I);
end

function I_thresh = ordered_thresh(I, M)

    % Implements the ordered threshold method. This is realized through an
    % if-else statement reduced to one line, i.e. (a > b)*c + (~(a > b))*d,
    % where we set c to 1 and d to 0. a is the value of the input image at
    % a certain position and b is the respective threshold value at that
    % position for the threshold mask.
    
    % param I:  input image to be thresholded
    % param M:  threshold mask
    % return:   binary image
    
    [height, width] = size(I);
    [height_m, width_m] = size(M);
    I_thresh = zeros(height, width);
    
    for x=0:width_m:(width - width_m - 1)
        for y=0:height_m:(height - height_m - 1)
            for x_m=1:width_m
                for y_m=1:height_m
                    % (a > b)*c + (~(a > b))*d
                    I_thresh(y + y_m, x + x_m) =  (I(y + y_m, x + x_m) > M(y_m, x_m)) * 1 + (~(I(y + y_m, x + x_m) > M(y_m, x_m))) * 0;
                end
            end
        end
    end
end

function I_quant = quantize(I, lower, upper)

    % Quantizes a grayscale image to a fixed range defined by an upper and
    % lower bound.
    
    % param I:      image to be quantized
    % param lower:  lower bound
    % param upper:  upper bound
    % return:       quantized image
    
    I_min = min(I, [], 'all');
    I_max = max(I, [], 'all');
    I_quant = uint8((((I - I_min)) .* (upper - lower)) ./ ((I_max - I_min) + lower));
end


function I = diffuse_pixel(I, x, y, error, c)

    % Diffuses one pixel of an image given the coordinates and a given
    % error value according to the error threshold method.
    %
    % param I:      working image
    % param x:      x-coordiante (centered pixel - offset)
    % param y:      y-coordinate (centered pixel - offset)
    % param error:  computed error for center pixel
    % param c:      coefficient for diffusing the error
    % return:       image with updated pixel at x, y
    
    I(y, x) = I(y, x) + error * c;
end


function I = error_diffusion(I, filter, t)

    % Implements the error diffusion method according to the lecture.
    %
    % param I:      input grayscale image
    % param filter: error diffusion filter of format y, x, coeff
    % param t:      threshold value
    % return:       dithered image
    
    [height, width] = size(I);
    for x=1:width
        for y=1:height
            old_pixel = I(y, x);
            
            % apply threshold on old pixel
            if (old_pixel >= t)
                new_pixel = 1.0;
            else 
                new_pixel = 0.0;
            end
            
            I(y, x) = new_pixel;
            % compute error
            error = old_pixel - new_pixel;
            
            for i=1:length(filter)
                % get offset coordinates
                y_i = y + filter(i, 1);
                x_i = x + filter(i, 2);
                
                % diffuse error
                if (1 <= y_i && y_i <= height) && (1 <= x_i && x_i <= width)
                    I = diffuse_pixel(I, x_i, y_i, error, filter(i, 3));
                end
            end
        end
    end
end