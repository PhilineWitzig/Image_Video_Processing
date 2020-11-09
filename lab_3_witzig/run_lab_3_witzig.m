% Lab 3 - Philine Witzig 07.11.2020

% getting path to project as it might vary from MATLAB user_path
path = matlab.desktop.editor.getActiveFilename;
path_split = strsplit(path, '/');
path_cur_folder = char(join(path_split(1:end - 1), '/'));

% choose the exercise to be executed
exercise = input("Enter the number of the exercise you want to execute: ");
I_lena = imread("Images/lena.png");
I_rice = imread("Images/rice.png");
I_road = imread("Images/road.png");


switch exercise
    case 1
    %% Exercise 1 - "Template Method"
    disp("Testing template methods...");
    sobel_g1 = 1/4 .* [1 0 -1; 2 0 -2; 1 0 -1];
    sobel_g2 = 1/4 .* [-1 -2 -1; 0 0 0; 1 2 1];
    
    prewitt_g1 = 1/3 .* [1 0 -1; 1 0 -1; 1 0 -1];
    prewitt_g2 = 1/3 .* [-1 -1 -1; 0 0 0; 1 1 1];
    
    roberts_g1 = [0 0 -1; 0 1 0; 0 0 0];
    roberts_g2 = [-1 0 0; 0 1 0; 0 0 0];
    
    % Tests on Lena image
    disp("Lena image, L1, Sobel")
    tic
    edge_img = edge_detect_template(double(I_lena), sobel_g1, sobel_g2, "L1", 25);
    toc
    figure('name', "Sobel, L1, lena.png");
    imshow(edge_img);
    
    % Testing robustness to noise
    I = imnoise(rescale(I_lena), 'gaussian', 5/255); % image has to be rescaled to [0,1], also sigma
    edge_img_noised = edge_detect_template(255 * I, sobel_g1, sobel_g2, "L1", 25);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=5): ", num2str(mse)));
    
    I = imnoise(rescale(I_lena), 'gaussian', 11/255);
    edge_img_noised = edge_detect_template(255 * I, sobel_g1, sobel_g2, "L1", 25);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=11): ", num2str(mse)));
    
    I = imnoise(rescale(I_lena), 'gaussian', 25/255);
    edge_img_noised = edge_detect_template(255 * I, sobel_g1, sobel_g2, "L1", 25);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=25): ", num2str(mse)));
    
    
    disp("Lena image, L2, Sobel")
    tic
    edge_img = edge_detect_template(double(I_lena), sobel_g1, sobel_g2, "L2", 25);
    toc
    figure('name', "Sobel, L2, lena.png");
    imshow(edge_img);
    
    
    I = imnoise(rescale(I_lena), 'gaussian', 5/255);
    edge_img_noised = edge_detect_template(255*I, sobel_g1, sobel_g2, "L2", 25);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=5): ", num2str(mse)));
    
    I = imnoise(rescale(I_lena), 'gaussian', 11/255);
    edge_img_noised = edge_detect_template(255*I, sobel_g1, sobel_g2, "L2", 25);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=11): ", num2str(mse)));
    
    I = imnoise(rescale(I_lena), 'gaussian', 25/255);
    edge_img_noised = edge_detect_template(255*I, sobel_g1, sobel_g2, "L2", 25);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=25): ", num2str(mse)));
    
    disp("Lena image, L2, Prewitt")
    tic
    edge_img = edge_detect_template(double(I_lena), prewitt_g1, prewitt_g2, "L2", 25);
    toc
    figure('name', "Prewitt, L2, lena.png");
    imshow(edge_img);
    
    I = imnoise(rescale(I_lena), 'gaussian', 5/255);
    edge_img_noised = edge_detect_template(255*I, prewitt_g1, prewitt_g2, "L2", 25);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=5): ", num2str(mse)));
    
    I = imnoise(rescale(I_lena), 'gaussian', 11/255);
    edge_img_noised = edge_detect_template(255*I, prewitt_g1, prewitt_g2, "L2", 25);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=11): ", num2str(mse)));
    
    I = imnoise(rescale(I_lena), 'gaussian', 25/255);
    edge_img_noised = edge_detect_template(255*I, prewitt_g1, prewitt_g2, "L2", 25);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=25): ", num2str(mse)));
    
    disp("Lena image, L2, Roberts")
    tic
    edge_img = edge_detect_template(double(I_lena), roberts_g1, roberts_g2, "L2", 25);
    toc
    figure('name', "Roberts, L2, lena.png");
    imshow(edge_img);
    
    I = imnoise(rescale(I_lena), 'gaussian', 5/255);
    edge_img_noised = edge_detect_template(255*I, roberts_g1, roberts_g2, "L2", 25);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=5): ", num2str(mse)));
    
    I = imnoise(rescale(I_lena), 'gaussian', 11/255);
    edge_img_noised = edge_detect_template(255*I, roberts_g1, roberts_g2, "L2", 25);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=11): ", num2str(mse)));
    
    I = imnoise(rescale(I_lena), 'gaussian', 25/255);
    edge_img_noised = edge_detect_template(255*I, roberts_g1, roberts_g2, "L2", 25);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=25): ", num2str(mse)));
    
    % Test on Rice image, using Sobel with L2 only
    disp("Rice image, L2, Sobel")
    tic
    edge_img = edge_detect_template(double(I_rice), sobel_g1, sobel_g2, "L2", 32);
    toc
    figure('name', "Sobel, L2, rice.png");
    imshow(edge_img);
    
    % Test on Road image, using Sobel with L2 only
    disp("Road image, L2, Sobel")
    tic
    edge_img = edge_detect_template(double(I_road), sobel_g1, sobel_g2, "L2", 30);
    toc
    figure('name', "Sobel, L2, road.png");
    imshow(edge_img);
    

    case 2
    %% Exercise 2 - "Compass Operator"
    disp("Testing compass operator...");
    
    kirsch_mask = [-3 -3 5; -3 0 5; -3 -3 5];
    
    % Test on Lena image
    disp("Lena image")
    tic
    edge_img = edge_detect_compass(double(I_lena), kirsch_mask, 25);
    toc
    figure('name', "Compass, lena");
    imshow(edge_img);
    
    % Testing noise 
    I = imnoise(rescale(I_lena), 'gaussian', 5/255);
    edge_img_noised = edge_detect_compass(255*I, kirsch_mask, 25);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=5): ", num2str(mse)));
    
    I = imnoise(rescale(I_lena), 'gaussian', 11/255);
    edge_img_noised = edge_detect_compass(255*I, kirsch_mask, 25);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=11): ", num2str(mse)));
    
    I = imnoise(rescale(I_lena), 'gaussian', 25/255);
    edge_img_noised = edge_detect_compass(255*I, kirsch_mask, 25);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=25): ", num2str(mse)));
    
    % Test on rice image
    disp("Rice image")
    tic
    edge_img = edge_detect_compass(double(I_rice), kirsch_mask, 32);
    toc
    figure('name', "Compass, rice");
    imshow(edge_img);
    
    %Test on road image
    disp("Road image")
    tic
    edge_img = edge_detect_compass(double(I_road), kirsch_mask, 25);
    toc
    figure('name', "Compass, road");
    imshow(edge_img);
    
    case 3
    %% Exercise 3 - "Laplace Operator"
    
    disp("LoG Operator")
    
    % Tests on Lena image, experimenting with sigma
    disp("Lena, sigma=2")
    tic
    edge_img = LoG(double(I_lena), 2);
    toc
    figure('name', "LoG, sigma=2");
    imshow(edge_img);
    
    disp("Lena, sigma=3")
    tic
    edge_img = LoG(double(I_lena), 3);
    toc
    figure('name', "LoG, sigma=3");
    imshow(edge_img);
    
    % Testing robustness to noise
    I = imnoise(rescale(I_lena), 'gaussian', 5/255);
    edge_img_noised = LoG(255*I, 3);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=5): ", num2str(mse)));
    
    I = imnoise(rescale(I_lena), 'gaussian', 11/255);
    edge_img_noised = LoG(255*I, 3);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=11): ", num2str(mse)));
    
    I = imnoise(rescale(I_lena), 'gaussian', 25/255);
    edge_img_noised = LoG(255*I, 3);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=25): ", num2str(mse)));
    
    disp("Lena, sigma=5")
    tic
    edge_img = LoG(double(I_lena), 5);
    toc
    figure('name', "LoG, sigma=5");
    imshow(edge_img);
    
    % Test on rice image
    disp("Rice, sigma=3")
    tic
    edge_img = LoG(double(I_rice), 3);
    toc
    figure('name', "LoG, sigma=3");
    imshow(edge_img);
    
    % Test on road image
    disp("Road, sigma=3")
    tic
    edge_img = LoG(double(I_road), 3);
    toc
    figure('name', "LoG, sigma=3");
    imshow(edge_img);
    
    
    case 4
    %% Exercise 4 - "Frei-Chen Method"
    
    disp("Frei-chen method");
    
    % Test on Lena image
    disp("Lena, thresh=0.97")
    tic
    edge_img = frei_chen(double(I_lena),0.97);
    toc
    figure('name', "Frei-Chen");
    imshow(edge_img);
    
    % Testing robustness to noise
    I = imnoise(rescale(I_lena), 'gaussian', 5/255);
    edge_img_noised = frei_chen(255*I,0.97);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=5): ", num2str(mse)));
    
    I = imnoise(rescale(I_lena), 'gaussian', 11/255);
    edge_img_noised = frei_chen(255*I,0.97);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=11): ", num2str(mse)));
    
    I = imnoise(rescale(I_lena), 'gaussian', 25);
    edge_img_noised = frei_chen(255*I,0.97);
    mse = immse(rescale(edge_img), double(edge_img_noised));
    disp(strcat("MSE noised image (sigma=25): ", num2str(mse)));
    
    % Test on rice image
    disp("Rice, thresh=0.98")
    tic
    edge_img = frei_chen(double(I_rice),0.98);
    toc
    figure('name', "Frei-Chen");
    imshow(edge_img);
    
    % Test on road image
    disp("Road, thresh=0.98")
    tic
    edge_img = frei_chen(double(I_road),0.98);
    toc
    figure('name', "Frei-Chen");
    imshow(edge_img);
    
end 


function edges = edge_detect_template(I, g1, g2, norm, thresh)
    % This function applies edge detection to an input image using the
    % template based method.
    
    % :param I:         input image in which we want to detect edges
    % :param g1:        template used for detecting vertical edges
    % :param g2:        template used for detecting horizontal edges
    % :param norm:      norm, either "L1" or "L2" used for approximating
    %                   the gradient
    % :param thresh:    threshold value in range (0, 1) for thresholding
    %                   filtered image
    % :return:          binary image, 1 for edge pixel, 0 for non-edge
    %                   pixel
    
    % filter image horizontally and vertically
    edges_vertical = conv2(I, g1, 'same');
    edges_horizontal = conv2(I, g2, 'same');
   
    % compute magnitude using respective norm
    if (norm == "L1" )
        magnitude = abs(edges_horizontal) + abs(edges_vertical);
       
    elseif (norm == "L2")
        magnitude = sqrt(edges_horizontal .^2 + edges_vertical .^2);
    else
        disp("Invalid norm for approximation.")
    end
    
    % extract edge pixels through thresholding
    edges = imbinarize(magnitude, thresh);
    
end

function edges = edge_detect_compass(I, mask, thresh)
    % Performs edge detection using the compass method, approximating the 
    % gradient in multiple different directions. The gradient magnitude is 
    % taken as the maximum absolute value amongst all directions.
    % 
    % :param I:         input image
    % :param mask:      reference mask which is going to be rotated for the
    %                   different directions
    % :param thresh:    threshold for extracting edge pixels
    
    % :return:          binary image, 1 for edge pixel, 0 for non-edge
    %                   pixel
    
    [height, width] = size(I);
    cur_max = zeros(height, width);
    
    for c=1:8
        % clockwise rotation of reference filter
        cur_mask = imrotate(mask, 45 * c, 'crop'); 
        cur_edges = abs(conv2(I, cur_mask, 'same'));
        % keep track of elementwise maximum pixel value
        cur_max = max(cur_max, cur_edges); 
    end
    
    % scale to [0, 255] for easy threshold comparison
    edges = 255 * rescale(cur_max);
    edges = imbinarize(edges, thresh);

end

function edges = LoG(I, sigma)
    % Performs edge detection using the Laplcaian of Gaussians. 
    % :param I:     input image
    % :param sigma: standard deviation
    
    % :return:      binary image, 1 for edge pixel, 0 for non-edge pixel
    
    [height, width] = size(I);
    edges = zeros(height, width);
    
    % LoG filter size
    hsize = 2 * ceil(3 * sigma) + 1;
    % getting the filter
    h = fspecial('log', hsize, sigma);
    % applying it through convolution
    I_filtered = conv2(I, h, 'same');
    
    % search for zero crossing by thresholding at 0
    % computing boundary pixels
    boundary_pxls = bwboundaries(imbinarize(I_filtered, 0));
    
    % color boundary pixels in white
    for k = 1 : length(boundary_pxls)
        cur_boundary = boundary_pxls{k};  % get k'th boundary
        edges(sub2ind(size(edges),cur_boundary(:,1),cur_boundary(:,2))) = 1;
    end
        
end

function edges = frei_chen(I, thresh)
    % Performs edge detection using the Frei-Chen method, i.e. creates a 9D
    % representation of the input image and then computes the angle of each
    % pixel between the 9D vector and its projection to the two edge
    % dimensions.
    
    % :param I:     input image
    % :param thesh: threshhold on angle values for separating edge from
    %               non-edge fixels
    % :return:      binary image, 1 for edge pixel, 0 for non-edge pixel
    [height, width] = size(I);
    f0 = 1 / sqrt(2) .* [1 sqrt(2) 1; 0 0 0; -1 -sqrt(2) -1];
    f1 = imrotate(f0, 90, 'crop');
    f2 = imrotate(f0, 135, 'crop');
    f3 = imrotate(f0, 215, 'crop');
    f4 = 1/2 .* [0 1 0; -1 0 -1; 0 1 0];
    f5 = imrotate(f4, 134, 'crop');
    f6 = 1/6 .* [1 -2 1; -2 4 -2; 1 -2 1];
    f7 = imrotate(f6, 45, 'crop');
    f8 = 1/3 .* ones(3, 3);
    
    I_9D = zeros(9, height, width);
  
    % disp(size(conv2(I, f0, 'same')));
    I_9D(1, :, :) = conv2(I, f0, 'same');
    I_9D(2, :, :) = conv2(I, f1, 'same'); 
    I_9D(3, :, :) = conv2(I, f2, 'same');
    I_9D(4, :, :) = conv2(I, f3, 'same');
    I_9D(5, :, :) = conv2(I, f4, 'same');
    I_9D(6, :, :) = conv2(I, f5, 'same');
    I_9D(7, :, :) = conv2(I, f6, 'same');
    I_9D(8, :, :) = conv2(I, f7, 'same');
    I_9D(9, :, :) = conv2(I, f8, 'same');
    
    m = zeros(height, width);
    s = zeros(height, width);
 
    for n=1:9
        cur_D = squeeze(I_9D(n, :, :));

        if (n <=2)
            m = m + cur_D.^2;
        end
        
        s = s + cur_D.^2;

    end
  
    angles = cos(sqrt(m ./ s));
    edges = imbinarize(angles, thresh);
    
end