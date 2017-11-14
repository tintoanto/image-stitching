function part2()
    %%
    %% load images and match files for the first example
    %%
    img_root = '../data/part2/';
    I1 = imread(strcat(img_root,'house1.jpg'));
    I2 = imread(strcat(img_root,'house2.jpg'));
    matches = load(strcat(img_root,'house_matches.txt')); 
    % this is a N x 4 file where the first two numbers of each row
    % are coordinates of corners in the first image and the last two
    % are coordinates of corresponding corners in the second image: 
    % matches(i,1:2) is a point in the first image
    % matches(i,3:4) is a corresponding point in the second image

    N = size(matches,1);

    %%
    %% display two images side-by-side with matches
    %% this code is to help you visualize the matches, you don't need
    %% to use it to produce the results for the assignment
    %%
%     imshow([I1 I2]); hold on;
%     plot(matches(:,1), matches(:,2), '+r');
%     plot(matches(:,3)+size(I1,2), matches(:,4), '+r');
%     line([matches(:,1) matches(:,3) + size(I1,2)]', matches(:,[2 4])', 'Color', 'r');
    % pause;

    %%
    %% display second image with epipolar lines reprojected 
    %% from the first image
    %%

    % first, fit fundamental matrix to the matches
    RANSAC      =   1;
    normalised  =   1;
    if RANSAC
        [ F, matches ] = ransac(I1,I2,normalised);
        N = size(matches,1);
    else
        F = fit_fundamental(matches,normalised); % this is a function that you should write
    end
    L = (F * [matches(:,1:2) ones(N,1)]')'; % transform points from 
    % the first image to get epipolar lines in the second image

    % find points on epipolar lines L closest to matches(:,3:4)
    L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
    pt_line_dist = sum(L .* [matches(:,3:4) ones(N,1)],2);
    closest_pt = matches(:,3:4) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);

    
    %% MSD pt line image 1
    pt_line_dist1_sq = pt_line_dist.^2;
    sum_pt_line_dist1_sq = sum(pt_line_dist1_sq);
    pt_line_dist1_sq_mean = sum_pt_line_dist1_sq/size(pt_line_dist1_sq,1);
    %%% MSD pt line image 2
    L2 = (F' * [matches(:,3:4) ones(N,1)]')'; % transform points from 
    L2 = L2 ./ repmat(sqrt(L2(:,1).^2 + L2(:,2).^2), 1, 3); % rescale the line
    pt_line_dist2 = sum(L2 .* [matches(:,1:2) ones(N,1)],2);
    closest_pt2 = matches(:,1:2) - L2(:,1:2) .* repmat(pt_line_dist2, 1, 2);

    pt_line_dist2_sq = pt_line_dist2.^2;
    sum_pt_line_dist2_sq = sum(pt_line_dist2_sq);
    pt_line_dist2_sq_mean = sum_pt_line_dist2_sq/size(pt_line_dist2_sq,1);
    
    pt3 = closest_pt2 - [L2(:,2) -L2(:,1)] * 10; % offset from the closest point is 10 pixels
    pt4 = closest_pt2 + [L2(:,2) -L2(:,1)] * 10;
    
    clf;
    imshow(I1); hold on;
    plot(matches(:,1), matches(:,2), '+r');
    line([matches(:,1) closest_pt2(:,1)]', [matches(:,2) closest_pt2(:,2)]', 'Color', 'r');
    line([pt3(:,1) pt4(:,1)]', [pt3(:,2) pt4(:,2)]', 'Color', 'g');
    %%
    
    fprintf("Image 1 points and epipolar of image 2 : " + pt_line_dist1_sq_mean +"\n")
    fprintf("Image 2 points and epipolar of image 1 : " + pt_line_dist2_sq_mean )

    
    % find endpoints of segment on epipolar line (for display purposes)
    pt1 = closest_pt - [L(:,2) -L(:,1)] * 10; % offset from the closest point is 10 pixels
    pt2 = closest_pt + [L(:,2) -L(:,1)] * 10;

    % display points and segments of corresponding epipolar lines
    clf;
    imshow(I2); hold on;
    plot(matches(:,3), matches(:,4), '+r');
    line([matches(:,3) closest_pt(:,1)]', [matches(:,4) closest_pt(:,2)]', 'Color', 'r');
    line([pt1(:,1) pt2(:,1)]', [pt1(:,2) pt2(:,2)]', 'Color', 'g');
    
    
    %% Load Matrix Camera and find Camera centers
    P1 = load(strcat(img_root,'house1_camera.txt'));
    P2 = load(strcat(img_root,'house2_camera.txt'));
    [~, ~, V] = svd(P1);
    cc1 = V(:,end);
    cc1 = cc1/cc1(4);
    cc1 = cc1(1:3);
    [~, ~, V] = svd(P2);
    cc2 = V(:,end);
    cc2 = cc2/cc2(4);
    cc2 = cc2(1:3);
    
    %% Triangulation
    for i = 1:size(matches,1)
        pt1 = matches(i,[1,2]);
        pt2 = matches(i,[3,4]);
        crossProductMat1 = [  0   -1  pt1(2); 1   0   -pt1(1); -pt1(2)  pt1(1)   0  ];
        crossProductMat2 = [  0   -1  pt2(2); 1   0   -pt2(1); -pt2(2)  pt2(1)   0  ];    
        Eqns = [ crossProductMat1*P1; crossProductMat2*P2 ];

        [~,~,V] = svd(Eqns);
        cal3D_homo = V(:,end)';
        cal3D_homo  = cal3D_homo/cal3D_homo(4);
        cal3D_cart(i,:) = cal3D_homo(1:3);

        projP1 = (P1 * cal3D_homo')';
        projP1 = projP1/projP1(3);
        projP1_cart(i,:) = projP1(1:2);
        projP2 = (P2 * cal3D_homo')';
        projP2 = projP2/projP2(3);
        projP2_cart(i,:) = projP2(1:2);
    end

    %% Residual Calculation
    res1_mat = pdist2(matches(:,1:2),projP1_cart);
    res1 = trace(res1_mat.^2);
    res1 = res1/size(res1_mat,1);
    fprintf("\n Residual Image 1 : "+ res1 );
    res2_mat = pdist2(matches(:,3:4),projP2_cart);
    res2 = trace(res2_mat.^2);
    res2 = res2/size(res1_mat,1);
    fprintf("\n Residual Image 2 : "+res2 );

    %% Plotting
    figure; axis equal;  hold on; 
    view(3);
    plot3(-cal3D_cart(:,1), cal3D_cart(:,2), cal3D_cart(:,3), '.r');
    plot3(-cc1(1), cc1(2), cc1(3),'*g');
    plot3(-cc2(1), cc2(2), cc2(3),'*b');
    grid on; xlabel('x'); ylabel('y'); zlabel('z'); axis equal;
    



    