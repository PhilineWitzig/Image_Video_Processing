% Lab 3 - Philine Witzig 03.11.2020

% getting path to project as it might vary from MATLAB user_path
path = matlab.desktop.editor.getActiveFilename;
path_split = strsplit(path, '/');
path_cur_folder = char(join(path_split(1:end - 1), '/'));

% choose the exercise to be executed
exercise = input("Enter the number of the exercise you want to execute: ");
I_lena = imread("Images/lena.png");


switch exercise
    case 1
    %% Exercise 1 - "Template Method"
    sobel_g1 = 1/4 .* [1 0 -1; 2 0 -2; 1 0 -1];
    sobel_g2 = 1/4 .* [-1 -2 -1; 0 0 0; 1 2 1];
    
    prewitt_g1 = 1/3 .* [1 0 -1; 1 0 -1; 1 0 -1];
    prewitt_g2 = 1/3 .* [-1 -1 -1; 0 0 0; 1 1 1];
    
    roberts_g1 = [0 0 -1; 0 1 0; 0 0 0];
    roberst_g2 = [-1 0 0; 0 1 0; 0 0 0];
    
    edge_img = edge_detect_template(double(I_lena), prewitt_g1, prewitt_g2, 25, "L1");
    figure('name', "Sobel");
    imshow(edge_img);
    

    case 2
    %% Exercise 2 - "Compass Operator"
    
    kirsch_mask = [-3 -3 5; -3 0 5; -3 -3 5];
    
    edge_img = edge_detect_compass(double(I_lena), kirsch_mask, 25);
    figure('name', "Compass");
    imshow(edge_img);
    
    case 3
    %% Exercise 3 - "Laplace Operator"
    
    edge_img = LoG(double(I_lena), 3);
    figure('name', "LoG");
    imshow(edge_img);

    case 4
    %% Exercise 4 - "Frei-Chen Method"
    edge_img = frei_chen(double(I_lena),0.99);
    figure('name', "Frei-Chen");
    imshow(edge_img);
    
end 


function edges = edge_detect_template(I, g1, g2, thresh, norm)
    % This function applies edge detection to an input image using the
    % template based method.
    
    % :param I:         input image in which we want to detect edges
    % :param template:  template used for filtering of form (gk, gl)
    % :param thresh:    threshold value in range (0, 1) for thresholding
    %                   filtered image
    % :param norm:      norm, either "L1" or "L2" used for approximating
    %                   the gradient
    % :return:          binary image, 1 for edge pixel, 0 for non-edge
    %                   pixel
    
    % filter padded image horizontally and vertically
    edges_vertical = conv2(I, g1, 'same');
    edges_horizontal = conv2(I, g2, 'same');
   
    if (norm == "L1" )
        magnitude = abs(edges_horizontal) + abs(edges_vertical);
       
    elseif (norm == "L2")
        magnitude = sqrt(edges_horizontal .^2 + edges_vertical .^2);
    else
        disp("Invalid norm for approximation.")
    end
    
    edges = imbinarize(magnitude, thresh);
    
end

function edges = edge_detect_compass(I, mask, thresh)
    [height, width] = size(I);
    cur_max = zeros(height, width);
    
    for theta=1:8
        cur_mask = imrotate(mask, 45*theta, 'crop'); % clockwise rotation but order does not matter
        cur_edges = abs(conv2(I, cur_mask, 'same'));
        cur_max = max(cur_max, cur_edges); 
    end
    
    edges = 255 * rescale(cur_max);
    edges = imbinarize(edges, thresh);

end

function edges = LoG(I, sigma)
    [height, width] = size(I);
    edges = zeros(height, width);
    hsize = 2 * ceil(3 * sigma) + 1;
    h = fspecial('log', hsize, sigma);
    I_filtered = conv2(I, h, 'same');
    boundary_pxls = bwboundaries(imbinarize(I_filtered, 0));
    for k = 1 : length(boundary_pxls)
        cur_boundary = boundary_pxls{k};  % Get k'th boundary
        edges(sub2ind(size(edges),cur_boundary(:,1),cur_boundary(:,2))) = 1;
    end
        
end

function edges = frei_chen(I, thresh)
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
    disp(size(I_9D(1, :, :)));
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
    disp(angles)
    edges = imbinarize(angles, thresh);
    
end