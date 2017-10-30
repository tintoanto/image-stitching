function main()
    IMatrix = preprocess();
    num = size(IMatrix,3);
    nsize = 5 ;                  % neighbourhood size
    

% Getting Corner coordinates
% Q2. Detect feature points in both images. You can use the Harris corner 
%     detector code harris.m that we provide or the blob detector you have 
%     developed as part of HW 2.

    r= {}; c= {};
    for i = 1:num
        [~, ri, ci] = harris(IMatrix(:,:,i), 2, 0.1, 2, 0);
        r(i) = {ri};
        c(i) = {ci};
    end
    
% Getting the neighborhoods
% Q3. Extract local neighborhoods around every keypoint in both images, 
%     and form descriptors simply by “flattening” the pixel values in each 
%     neighborhood to one-dimensional vectors. Experiment with different 
%     neighborhood sizes to see which one works the best. If you’re using your 
%     Laplacian detector, use the detected feature scales to define the 
%     neighborhood scales.
    n= {};
    for i = 1:num
        ni = neighborhoods(IMatrix(:,:,i), r{i}, c{i}, nsize);
        n(i)={ni};
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
    disp(size(d))

% Select the matches
% Q5. Select putative matches based on the matrix of pairwise descriptor 
%     distances obtained above. You can select all pairs whose descriptor 
%     distances are below a specified threshold, or select the top few hundred
%     descriptor pairs with the smallest pairwise distances.
    [M, dJ] = min(d,[],2);
    [B,oriDJ] = sort(M);
    
    



























