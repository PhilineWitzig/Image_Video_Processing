% Lab 1 - Philine Witzig and 06.10.2020

% getting path to project as it might vary from MATLAB user_path
path = matlab.desktop.editor.getActiveFilename;
path_split = strsplit(path, '/');
path_cur_folder = char(join(path_split(1:end - 1), '/'));

% choose the exercise to be executed
exercise = input("Enter the number of the exercise you want to execute: ");

switch exercise
    case 1
        %% Exercise 1 - "Images and Color Tables"
        disp("Ex. 1 ...")
        
        % load images
        tree_file = fullfile(path_cur_folder, 'Images', 'trees.tif');
        lena_file = fullfile(path_cur_folder, 'Images', 'lena.tif');
        [I_trees, map_trees] = imread(tree_file);
        I_lena = imread(lena_file);

        % display original images
        disp("Original images");

        imshow(I_trees, map_trees);
        title("Trees original");
        pause();

        imshow(I_lena);
        title("Lena original");
        pause();

        % display images in gray scale
        disp("Images in grayscale");
        map_trees_gray = rgb2gray(map_trees);
        I_lena_gray = rgb2gray(I_lena);

        imshow(I_trees, map_trees_gray);
        title("Trees grayscale");
        pause();

        imshow(I_lena_gray);
        title("Lena grayscale");
        pause();

        % invert images
        disp("Inverted images")
        map_trees_gray_inv = imcomplement(map_trees_gray);
        I_lena_gray_inv = imcomplement(I_lena_gray);

        imshow(I_trees, map_trees_gray_inv);
        title("Trees grayscale inverted");
        pause();

        imshow(I_lena_gray_inv);
        title("Lena grayscale inverted");
        pause();

        % gamma correction
        disp("Gamma corrected images");
        gamma_values = [0.5 0.75 1.25 1.5 1.75 2.0]; % gamma test candidates

        for i=1:length(gamma_values)
            map_trees_gamma = gamma_corr(map_trees, gamma_values(i));
            imshow(I_trees, map_trees_gamma);
            title(strcat("Trees gamma corrected, gamma = ", num2str(gamma_values(i))));
            pause();
        end

        % create chess board
        disp("Creating chess board.")
        I_chess = repmat([1 2; 2 1], 4, 4);
        cmap = [1 1 0; 0 0 1]; % [yellow, blue]
        map_chess = colormap(cmap);

        path_chess_indexed = fullfile(path_cur_folder, 'Images', 'chess_indexd.tif');
        path_chess_rgb = fullfile(path_cur_folder, 'Images', 'chess_rgb.tif');
        imwrite(I_chess, map_chess, path_chess_indexed);
        imwrite(ind2rgb(I_chess, map_chess), path_chess_rgb); % convert ind to rgb and save it

        
    case 2
        %% Exercise 2 - "Image Qunatization"
        disp('Ex. 2 ...')
        
        lena_y_file = fullfile(path_cur_folder, 'Images', 'lena-y.png');
        I_lena_y = imread(lena_y_file);
        levels = [128 64 32 16 8 4 2];

        for i=1:length(levels)
            bin_size = 256 / (levels(i) - 1);
            I_lena_y_quant = bin_size * floor(I_lena_y ./ (bin_size));
            imshow(I_lena_y_quant);
            title(strcat("Uniform quantization of lena with ", num2str(levels(i)), " levels."));
            pause();
            cur_path = fullfile(path_cur_folder, 'Images', strcat('lena-y_', num2str(levels(i)), '.png'));
            imwrite(I_lena_y_quant, cur_path);
        end
        
        
    case 3
        %% Exercise 3 - "Filtering"
        disp("Ex. 3 ...")
        I_gold_text = imread(fullfile(path_cur_folder, 'Images', 'gold-text.png'));
        
        % 5x5 filter
        F1 = [0.0357; 0.2411; 0.4464; 0.2411; 0.0357] * [0.0357 0.2411 0.4464 0.2411 0.0357];
        freqz2(F1);
        title("Frequency response of 5x5 Filter");
        pause();
        
        I_gold_text_filt = imfilter(I_gold_text, F1, 'conv');
        imshow(I_gold_text_filt, []);
        title("Test image convolved with 5x5 filter, edges trimmed");
        pause();
        
        % 3x3 filter
        F2 = (1 / 6) .* [-1 -4 -1; -4 26 -4; -1 -4 -1];
        freqz2(F2);
        title("Frequence response of 3x3 Filter");
        pause();
        I_gold_text_filt = imfilter(I_gold_text, F2, 'conv');
        imshow(I_gold_text_filt);
        title("Test image convolved with 3x3 filter, edges trimmed");
        pause();

        
    case 4
        %% Exercise 4 - "Correlation"
        disp("Ex. 4 ...")
        disp("Correlation in spatial domain")
        I_gold_text = double(imread(fullfile(path_cur_folder, 'Images', 'gold-text.png')));
        I_gold_text = I_gold_text - 128; % make image have zero nominal 4average
        template = double(imread(fullfile(path_cur_folder, 'Images', 'g-letter.png')));
        template = template - 128; % make template have zero nominal average
        phi = conv2(I_gold_text, fliplr(template));
        [height, width] = size(I_gold_text);
       
        % visualize max response
        mask = get_mask(I_gold_text, template, phi);
        imshow(mask);
        title("Ex. 4: max. response after correlation (spatial)");
        pause();

        % perform correlation in the frequency domain
        disp("Correlation in frequency domain");
        I_freq = fft2(I_gold_text);
        temp_freq = fft2(template);
        phi_freq = conv2(I_freq, fliplr(temp_freq));
       
        mask = get_mask(I_gold_text, template, ifft2(phi_freq));
        imshow(mask);
        title("Ex. 4: max. response after correlation (frequency)");
        pause();
        
        % add noise
        disp("Correlation in spatial domain with noise added to image");
        stdvs = [5 10 25 40 50];
        for i=1:length(stdvs)
            noise = uint8(255 * mat2gray(stdvs(i) .* randn(height, width))); % create noise with stdv
            noise = double(noise) - 128; % shift to range -128 to 127
            I_gold_text_n = double(I_gold_text) + noise;
            phi = conv2(I_gold_text_n, fliplr(template));
            mask = get_mask(I_gold_text_n, template, phi);
            imshow(mask);
            title(strcat("Ex. 4: max. response after correlation (spatial) width added noise (stdv = ", num2str(stdvs(i)), ')'));
            pause();
        end
        
        % normalized cross correlation
        disp("Investigating normalized cross correlation");
        phi = normxcorr2(template, I_gold_text);
        mask = get_mask(I_gold_text_n, template, phi);
        imshow(mask);
        title("Ex. 4: normalized cross correlation");

    case 5
        %% Exercise 5 - "Resampling"
        disp("Ex. 5 ...");
        I = imread(fullfile(path_cur_folder, 'Images', 'sub4.tif'));
        I_sub_2 = I(1:2:end, 1:2:end); % downsample by factor 2
        I_sub_4 = I(1:4:end, 1:4:end); % downsample by factor 4
        
        imshow(I);
        title("Ex. 5: original image");
        pause();
        imshow(I_sub_2);
        title("Ex. 5: downsampled, factor 2");
        pause();
        imshow(I_sub_4);
        title("Ex. 5: downsampled, factor 4");
        pause();

    case 6
        %% Exercise 6 - "Phase and Magnitude of th 2DFT"
        disp("Ex.6 ...");
        
        I_lena = imread(fullfile(path_cur_folder, 'Images', 'lena-y.png'));
        I_lena_FT = fft2(I_lena);
        
        % 2DIFT without imaginary part
        I_lena_real = ifft2(real(I_lena_FT)); % result is real
        imshow(I_lena_real, []);
        title("Ex.6: 2DIFT of lena without imag. part");
        pause();
        
        % 2DIFT without real part
        I_lena_imag = ifft2(imag(I_lena_FT)); % result is complex
        imshow(real(I_lena_imag), []);
        title("Ex.6: 2DIFT of lena without real part");
        pause();
        
        % 0-phase
        I_lena_FT_0_phase = abs(I_lena_FT); % cf. lecture
        imshow(log(real(ifft2(I_lena_FT_0_phase))), []);
        title("Ex.6: 2DIFT with 0 phase");
        pause();
        
        % magnitude of 1
        I_lena_FT_mag_1 = exp(1i*angle(I_lena_FT)); % cf. lecture
        imshow(real(ifft2(I_lena_FT_mag_1)), []);
        title("Ex.6: 2DIFT with magnitude 1");
        pause();
        
    case 7
        %% Exercise 7 - "Weber Law"
        disp("Ex. 7: ...")
        disp("Estimating your Weber constant.");
        disp("Press any button to switch from image display to console input to answer questions.");
        diff = true;
        lb = 10;
        l1 = 230;
        l2 = 250;
        l1_last = l1;
        l2_last = l2;
        incr = 10;
        
        while (diff == true)
            I = weber(l1, l2, lb);
            imshow(I);
            pause();
            resp = input("Did you notice a difference? [y/n]", 's');
            % user response is positive
            if (resp == 'y')
                l1_last = l1;
                l2_last = l2;
                % choose correct increment
                if (l1 + 20 < l2)
                    incr = 10;
                elseif (l1 + 10 < l2)
                    incr = 5;
                else
                    incr = 1;
                end
                l1 = l1 + incr;
                l2 = l2 - incr;
            % user response is negative
            elseif (resp == 'n')
                if (incr == 1)
                    diff = false;
                elseif (incr == 5)
                    incr = 1;
                    l1 = l1_last;
                    l2 = l2_last;
                elseif (incr == 10)
                    incr = 5;
                    l1 = l1_last;
                    l2 = l2_last;
                else
                    disp("Unhandled condition.")
                end
                
            end
        end
        
        disp(strcat("Weber constant is: ", num2str((l2 - l1) / l1), ". Background intensity value: ", num2str(lb)));
        

    otherwise
        disp("Invalid exercise number.")
        
end

%% Helper functions
function pause()
    % pausing image display until any key pressed
    try
        waitforbuttonpress;
    catch
        disp("Image closed.");
    end
end

function new_color_table = gamma_corr(color_table, gamma)
    % gamma correction
    new_color_table = color_table.^gamma; % elementwise exponential, gamma is the same for all channels
end

function mask = get_mask(I, template, phi)
    % visualize max response of correlation
    phi_max = max(phi, [], 'all');
    [y, x] = find(phi == phi_max);
    offset = [y - size(template, 1), x - size(template, 2)];
    x_begin = offset(2) + 1;
    x_end = offset(2) + size(template, 2);
    y_begin = offset(1) + 1;
    y_end = offset(1) + size(template, 1);
    mask = uint8(I + 128);
    mask(y_begin: y_end, x_begin: x_end) = uint8(template + 128);
end

function I = weber(l1, l2, lb)
    % produce weber test image
    I = lb .* ones(600,600);
    I(220:380, 220:300) = l1;
    I(220:380, 300:380) = l2;
    I = uint8(I);
end