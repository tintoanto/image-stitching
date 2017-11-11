function [ total_error, inlier_error, inliers, inlier_count ] = find_error(pm, pm_idx, h, threshold)
    total_error     = 0; 
    inlier_error    = 0;
    inliers         = [];
    inlier_count    = 0;
    
    im1x    = pm{1}(:,1);
    im1y    = pm{1}(:,2);
    x1      = im1x(pm_idx);
    y1      = im1y(pm_idx);

    im2x    = pm{2}(:,1);
    im2y    = pm{2}(:,2);
    x2      = im2x(pm_idx);
    y2      = im2y(pm_idx);
    
    for i = 1:size(pm_idx,2)
        p1 = [ x1(i); y1(i); 1];
        p2 = [ x2(i); y2(i); 1]';
        p2_predict =    h*p1;
        p2_predict = (p2_predict./p2_predict(3))';
        error = pdist2(p2(1:2),p2_predict(1:2));
        total_error = total_error + error;
        if error< threshold
            inlier_error    =   inlier_error + error;
            inliers         =   [ inliers; pm_idx(i) ];
            inlier_count    =   inlier_count + 1;
        end
    end
    
%   disp(inlier_count);
    
    