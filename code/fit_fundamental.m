function F = fit_fundamental(matches,normalised)
%     samples = randperm(size(matches,1));
%     matches = matches(samples,1:4);


    if normalised
        centroid = sum(matches(1:8,1:2)) / 8;
        diff = matches(1:8,1:2) - repmat(centroid,8,1);
        avgdist = sum(sqrt(diff(:,1).^2 + diff( :,2).^2)) / 8;
        scale = sqrt(2) / avgdist;
        T = diag([scale scale 1]) * [eye(3,2) [-centroid 1]'];
        newpts = matches(1:8,1:2);
        newpts(:,3) = 1;
        newpts = newpts * T';
        newpts = newpts ./ repmat( newpts(:,3), 1, 3);
        x1 = newpts(:,1);
        y1 = newpts(:,2);
        
        
        centroid2 = sum(matches(1:8,3:4)) / 8;
        diff2 = matches(1:8,3:4) - repmat(centroid2,8,1);
        avgdist2 = sum(sqrt(diff2(:,1).^2 + diff2( :,2).^2)) / 8;
        scale2 = sqrt(2) / avgdist2;
        T2 = diag([scale2 scale2 1]) * [eye(3,2) [-centroid2 1]'];
        newpts2 = matches(1:8,3:4);
        newpts2(:,3) = 1;
        newpts2 = newpts2 * T2';
        newpts2 = newpts2 ./ repmat( newpts2(:,3), 1, 3);
        x2 = newpts2(:,1);
        y2 = newpts2(:,2);
    else
        x1  = matches(1:8,1);
        y1  = matches(1:8,2);
        x2  = matches(1:8,3);
        y2  = matches(1:8,4);        
    end
    
    
    A = [ x2.*x1 x2.*y1 x2 y2.*x1 y2.*y1 y2 x1 y1 ones(8,1)];
    [U,D,V] = svd(A,0);
    f = V(:,end);
    F = reshape( f,[3 3])';
    
    [U,D,V] = svd(F);
    D(3,3) = 0;             
    D = D / D(1,1);
    F = U * D * V';
    
    if normalised
        F = T2' * F * T;
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    