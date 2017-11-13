function [ F, matches ] = ransac(ori_I1,ori_I2,normalised)
    gray_I1  = rgb2gray(ori_I1);
    I1       = im2double(gray_I1);
    gray_I2  = rgb2gray(ori_I2);
    I2       = im2double(gray_I2);
    
    % Harris Params
    sigma           = 1;
    threshold       = 0.02;
    radius          = 1;
    display_harris  = 0;
    
    % Putatuive Matches
    nsize = 3;  
    pm_threshold = 0.01;
%     display_putative_matches = 1;

%     % Ransac
    ransac_count = 10000;
    ransac_threshold= 1000;
%     display_ransac_matches = 1;
% 
% 
%    
%% 

    [~, r1, c1]  = harris(I1, sigma, threshold, radius, display_harris);   
    [~, r2, c2]  = harris(I2, sigma, threshold, radius, display_harris);   
    
    [n1, ~, ~] = neighborhoods(I1, r1, c1, nsize);
    [n2, ~, ~] = neighborhoods(I2, r2, c2, nsize);
    
    d = dist2(n1,n2);   
    disp(size(d));
    [dist, pj] = min(d,[],2);
    pi = (1:size(d,1))';
    pi = pi(dist<pm_threshold);
    pj = pj(dist<pm_threshold);
    
    pmatches = [c1(pi) r1(pi) c2(pj) r2(pj)];
 
    
%%  RANSAC 
    N = size(pmatches,1);
    max_inlier_count = 0;
    best_F =[];
    best_inliers_idx = [];
    pmatches_idx = (1:size(pmatches,1))';
    
    for i = 1:ransac_count
        F = fit_fundamental(pmatches,normalised);
        
        L = (F * [pmatches(:,1:2) ones(N,1)]')'; % transform points from 
        % the first image to get epipolar lines in the second image

        % find points on epipolar lines L closest to matches(:,3:4)
        L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
        pt_line_dist = sum(L .* [pmatches(:,3:4) ones(N,1)],2);
%         closest_pt = matches(:,3:4) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);
        pt_line_dist1_sq = pt_line_dist.^2;
        sum_pt_line_dist1_sq = sum(pt_line_dist1_sq);
        pt_line_dist1_sq_mean = sum_pt_line_dist1_sq/size(pt_line_dist1_sq,1);
        root_pt_line_dist1_sq_mean = sqrt(pt_line_dist1_sq_mean);
        
        %%% MSD pt line image 2
        L2 = (F' * [pmatches(:,3:4) ones(N,1)]')'; % transform points from 
        L2 = L2 ./ repmat(sqrt(L2(:,1).^2 + L2(:,2).^2), 1, 3); % rescale the line
        pt_line_dist2 = sum(L2 .* [pmatches(:,1:2) ones(N,1)],2);
%          closest_pt2 = matches(:,1:2) - L2(:,1:2) .* repmat(pt_line_dist2, 1, 2);
        pt_line_dist2_sq = pt_line_dist2.^2;
        sum_pt_line_dist2_sq = sum(pt_line_dist2_sq);
        pt_line_dist2_sq_mean = sum_pt_line_dist2_sq/size(pt_line_dist2_sq,1);
        root_pt_line_dist2_sq_mean = sqrt(pt_line_dist2_sq_mean);

        inliers_bool = (pt_line_dist1_sq<ransac_threshold)&(pt_line_dist2_sq<ransac_threshold);
        inliers_count = sum(inliers_bool);
        inliers_idx = pmatches_idx((pt_line_dist1_sq<ransac_threshold)&(pt_line_dist2_sq<ransac_threshold));
        if(inliers_count > max_inlier_count)
            max_inlier_count=inliers_count;
            best_F = F;
            best_inliers_idx = inliers_idx;
        end
        
        error = root_pt_line_dist1_sq_mean + root_pt_line_dist2_sq_mean;
        
    end
    
    fprintf("Inliers : " + max_inlier_count + "\n");
    matches = pmatches( best_inliers_idx, 1:4 );
    F = best_F;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    