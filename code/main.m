function main()
    [IMatrix img1 img2]= preprocess(4);           % selecting folder and setting all images to IMatrix
                                       % IMatrix(:,:,idx) denotes 'idx'th image matrix
    num = size(IMatrix,3);             % number of images
    nsize = 13;                         % neighbourhood size

% Image(1)
%     % Harris Params
%     sigma       = 1;
%     threshold   = 0.075;
%     radius      = 1;
%     display     = 0;
%     
%     % Putatuive Matches
%     pm_count = 400;
% 
%     % Ransac
%     ransac_count = 400;
%     ransac_threshold= 0.65;

    % Harris Params
    sigma       = 1;
    threshold   = 0.15;
    radius      = 1;
    display_harris     = 0;
    
    % Putatuive Matches
    pm_count = 100;
    display_putative_matches = 1;

    % Ransac
    ransac_count = 100;
    ransac_threshold= 1;
    display_ransac_matches = 1;

    % Stitching
    display_final =1;
%% Getting Corner coordinates
%  Q2. Detect feature points in both images. You can use the Harris corner 
%     detector code harris.m that we provide or the blob detector you have 
%     developed as part of HW 2.

    r= {}; c= {};
    for idx = 1:num
        [~, ridx, cidx]  = harris(IMatrix(:,:,idx), sigma, threshold, radius, display_harris);   
        I(idx) = {ridx};                % I(idx) is vector denotes the I coordinates features of matrix idx
        J(idx) = {cidx};                % J(idx) is vector  denotes the J coordinates features of matrix idx
    end
    
%% Getting the neighborhoods
%  Q3. Extract local neighborhoods around every keypoint in both images, 
%     and form descriptors simply by “flattening” the pixel values in each 
%     neighborhood to one-dimensional vectors. Experiment with different 
%     neighborhood sizes to see which one works the best. If you’re using your 
%     Laplacian detector, use the detected feature scales to define the 
%     neighborhood scales.
    n= {};newI={};newJ={};
    for i = 1:num
        [ni newIi newJi] = neighborhoods(IMatrix(:,:,i), I{i}, J{i}, nsize);
        newI(i) = {newIi};              % newI(i) is reduced I(i)
        newJ(i) = {newJi};              % newJ(i) is reduced J(i)
        n(i)    = {ni};                 %  n(i) denotes the i'th features neighbourhood
    end
    
%% Get distances
%  Q4. Compute distances between every descriptor in one image and every 
%     descriptor in the other image. You can use dist2.m that we provide for fast 
%     computation of Euclidean distance. Alternatively, experiment with computing 
%     normalized correlation, or Euclidean distance after normalizing all 
%     descriptors to have zero mean and unit standard deviation. Optionally, feel 
%     free to experiment with SIFT descriptors. The script find sift.m that we 
%     provide contains some basic code for computing SIFT descriptors of
%     circular regions, such as the ones returned by the detector from HW 2.
    d = dist2(n{1},n{2});   
    disp(size(d));

%% Select the matches
%  Q5. Select putative matches based on the matrix of pairwise descriptor 
%     distances obtained above. You can select all pairs whose descriptor 
%     distances are below a specified threshold, or select the top few hundred
%     descriptor pairs with the smallest pairwise distances.
    
    [M, fulldJ] = min(d,[],2);
    [B, fulldI] = sort(M);
    dI          = fulldI(1:pm_count);
    dJ          = fulldJ(dI);
    
    
    pm          = {};                                              % pm = putative matches
    pm(1)       = { [ newJ{1}(dI) newI{1}(dI) ] };
    pm(2)       = { [ newJ{2}(dJ) newI{2}(dJ) ] };
    x_offset    = size(IMatrix(:,:,1),2);
    y_offset    = size(IMatrix(:,:,1),1);
    
    % display
    if display_putative_matches
        imshow( [   IMatrix(:,:,1) IMatrix(:,:,2)   ] ); hold on;
        scatter(    pm{1}(:,1),             pm{1}(:,2))
        scatter(    pm{2}(:,1)+x_offset,    pm{2}(:,2))
        line(   [   pm{1}(:,1)  pm{2}(:,1)+x_offset]',  ...
                [   pm{1}(:,2)  pm{2}(:,2)]',           ...
                        'Color', 'r');
    end
    
    
%% Run RANSAC
%  Q6. Run RANSAC to estimate a homography mapping one image onto the other. Report the number of
%     inliers and the average residual for the inliers (squared distance between the point coordinates in one
%     image and the transformed coordinates of the matching point in the other image). Also, display the
%     locations of inlier matches in both images.

    max_inlier_count = 0;
    best_h =[];
    best_inliers = [];
    best_inlier_error = 0;
    best_pm_random_4_idx =[];
    
    for i = 1:ransac_count
        pm_random_idx       =   randperm(pm_count);
        pm_random_4_idx       =   pm_random_idx(1:4);

        h = get_homography_vector(pm, pm_random_4_idx);

        pm_random_rest_idx       =   pm_random_idx(5:pm_count);
        
        [total_error inlier_error inliers inlier_count] = find_error(pm, pm_random_rest_idx, h, ransac_threshold);
        
        if inlier_count > max_inlier_count
            max_inlier_count = inlier_count;
            best_h = h;
            best_inliers = inliers;
            best_inlier_error = inlier_error;
            best_pm_random_4_idx=pm_random_4_idx;
        end
    end
    
    fprintf("Inliers : " + max_inlier_count +...
        "\nAvg Residual = " +  best_inlier_error/max_inlier_count  +...
        "\nHomography Matrix : \n")
    disp( best_h);

    rm  = {};
    
    px1_temp    = pm{1}(:,1);
    px1         = [px1_temp(best_pm_random_4_idx); px1_temp(best_inliers)];
    py1_temp    = pm{1}(:,2);
    py1         = [py1_temp(best_pm_random_4_idx); py1_temp(best_inliers)];
    
    px2_temp    = pm{2}(:,1);
    px2         = [px2_temp(best_pm_random_4_idx); px2_temp(best_inliers)];
    py2_temp    = pm{2}(:,2);
    py2         = [py2_temp(best_pm_random_4_idx); py2_temp(best_inliers)];
    
    rm(1)   =   {   [   px1     py1    ]   };
    rm(2)   =   {   [   px2     py2    ]   };
    
    % display
    if display_ransac_matches
        imshow( [   IMatrix(:,:,1) IMatrix(:,:,2)   ] ); hold on;
        line(   [   px1  px2 + x_offset     ]',  ...
                [   py1  py2                ]',           ...
                        'Color', 'r');    
        line(   [   rm{1}(:,1)  rm{2}(:,1)+x_offset]',  ...
                [   rm{1}(:,2)  rm{2}(:,2)]',           ...
                        'Color', 'r');
    end

%% Transform    
%  Q7. Warp one image onto the other using the estimated transformation. To do this, you will need to learn
%     about maketform and imtransform functions.
    T = maketform('projective',pinv(best_h'));
    [im2t,xdataim2t,ydataim2t] = imtransform(IMatrix(:,:,2),T);
    
%% Image stitching    
%  Q8. Create a new image big enough to hold the panorama and composite the two images into it. You can
%     composite by simply averaging the pixel values where the two images overlap. Your result should look
%     something like this (but hopefully with a more precise alignment):
    
    xdataout=[min(1,xdataim2t(1)) max(size(IMatrix(:,:,1),2),xdataim2t(2))];
    ydataout=[min(1,ydataim2t(1)) max(size(IMatrix(:,:,1),1),ydataim2t(2))];

    % let's transform both images with the computed xdata and ydata
    im2t=imtransform(img2,T,'XData',xdataout,'YData',ydataout);
    im1t=imtransform(img1,maketform('affine',eye(3)),'XData',xdataout,'YData',ydataout);

    
    
    if display_final
        ims = im1t/2+im2t/2;
        idx = ( 1 : size(im1t,1) * size(im1t,2) * size(im1t,3) );
        idx = idx( im1t==0 | im2t==0 );
        ims(idx) = ims(idx).*2;
        figure, imshow(ims)
    end
%%
    
    



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    