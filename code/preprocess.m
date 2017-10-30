function filesGrayScaleMatrix = preprocess(index)
%  Q1. Load both images, convert to double and to grayscale.
    if nargin < 1
        index = 4;
    end
    
    dataroot = '../data/part1';
    folderDescriptorsAll = dir(dataroot);
    folderDescriptors = folderDescriptorsAll(3:length(folderDescriptorsAll),:);
    
    folderName = folderDescriptors(index).name;
    folderPath = strcat(dataroot,'/', folderName);
    
    fileDescriptorsAll = dir(folderPath);
    fileDescriptors = fileDescriptorsAll(3:length(fileDescriptorsAll),:);
    
    filesGrayScaleMatrix = [];
    for i=1:length(fileDescriptors)
        fileName    = fileDescriptors(i).name;
        filePath    = strcat(folderPath,'/',fileName);

        ori_I   = imread(filePath);
        gray_I  = rgb2gray(ori_I);
        I       = im2double(gray_I);

        filesGrayScaleMatrix = cat (3,filesGrayScaleMatrix,I);
    end
    
    
