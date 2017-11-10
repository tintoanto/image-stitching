function main()
    IMatrix = preprocess(1);
    num = size(IMatrix,3);
    nsize = 9;                         % neighbourhood size
    

% Getting Corner coordinates
% Q2. Detect feature points in both images. You can use the Harris corner 
%     detector code harris.m that we provide or the blob detector you have 
%     developed as part of HW 2.

    r= {}; c= {};
    for idx = 1:num
        % Usage:  [cim, r, c] = harris(im,               sigma, thresh, radius, disp)
        [~, ridx, cidx]       = harris(IMatrix(:,:,idx), 1,     0.075,    1,      0);   % IMatrix(:,:,idx) denotes 'idx'th image matrix
        I(idx) = {ridx};                % I(idx) is vector denotes the I coordinates features of matrix idx
        J(idx) = {cidx};                % J(idx) is vector  denotes the J coordinates features of matrix idx
    end
    
% Getting the neighborhoods
% Q3. Extract local neighborhoods around every keypoint in both images, 
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
    
% Get distances
% Q4. Compute distances between every descriptor in one image and every 
%     descriptor in the other image. You can use dist2.m that we provide for fast 
%     computation of Euclidean distance. Alternatively, experiment with computing 
%     normalized correlation, or Euclidean distance after normalizing all 
%     descriptors to have zero mean and unit standard deviation. Optionally, feel 
%     free to experiment with SIFT descriptors. The script find sift.m that we 
%     provide contains some basic code for computing SIFT descriptors of
%     circular regions, such as the ones returned by the detector from HW 2.
    d = dist2(n{1},n{2});   
     disp(size(d));

% Select the matches
% Q5. Select putative matches based on the matrix of pairwise descriptor 
%     distances obtained above. You can select all pairs whose descriptor 
%     distances are below a specified threshold, or select the top few hundred
%     descriptor pairs with the smallest pairwise distances.
    pm_count = 400;
    [M, fulldJ] = min(d,[],2);
    [B, fulldI] = sort(M);
    dI          = fulldI(1:pm_count);
    dJ          = fulldJ(dI);
    
    
    pm          = {};                                              % pm = putative matches
    pm(1)       = { [ newJ{1}(dI) newI{1}(dI) ] };
    pm(2)       = { [ newJ{2}(dJ) newI{2}(dJ) ] };
    x_offset    = size(IMatrix(:,:,1),2);
    
    % display
%     imshow( [   IMatrix(:,:,1) IMatrix(:,:,2)   ] ); hold on;
%     scatter(    pm{1}(:,1),             pm{1}(:,2))
%     scatter(    pm{2}(:,1)+x_offset,    pm{2}(:,2))
%     line(   [   pm{1}(:,1)  pm{2}(:,1)+x_offset]',  ...
%             [   pm{1}(:,2)  pm{2}(:,2)]',           ...
%                     'Color', 'r');
    
    
    
% Run RANSAC
% Q6. Run RANSAC to estimate a homography mapping one image onto the other. Report the number of
%     inliers and the average residual for the inliers (squared distance between the point coordinates in one
%     image and the transformed coordinates of the matching point in the other image). Also, display the
%     locations of inlier matches in both images.

    max_inlier_count = 0;
    best_h =[];
    best_inliers = [];
    best_inlier_error = 0;
    
    for i = 1:400
        pm_random_idx       =   randperm(pm_count);
        pm_random_4_idx       =   pm_random_idx(1:4);

        h = get_homography_vector(pm, pm_random_4_idx);

        pm_random_rest_idx       =   pm_random_idx(5:pm_count);
        threshold= 0.65;
        [total_error inlier_error inliers inlier_count] = find_error(pm, pm_random_rest_idx, h, threshold);
        
        if inlier_count > max_inlier_count
            max_inlier_count = inlier_count;
            best_h = h;
            best_inliers = inliers;
            best_inlier_error = inlier_error;
        end
    end
    
    disp(max_inlier_count);
    disp(best_h);



    rm  = {};
    
    px1_temp    = pm{1}(:,1);
    px1         = [px1_temp(pm_random_4_idx); px1_temp(best_inliers)];
    py1_temp    = pm{1}(:,2);
    py1         = [py1_temp(pm_random_4_idx); py1_temp(best_inliers)];
    
    px2_temp    = pm{2}(:,1);
    px2         = [px2_temp(pm_random_4_idx); px2_temp(best_inliers)];
    py2_temp    = pm{2}(:,2);
    py2         = [py2_temp(pm_random_4_idx); py2_temp(best_inliers)];
    
    

    
    % display
    imshow( [   IMatrix(:,:,1) IMatrix(:,:,2)   ] ); hold on;
    line(   [   px1  px2 + x_offset     ]',  ...
            [   py1  py2                ]',           ...
                    'Color', 'r');    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    