function [neighborhood, newr, newc] = neighborhoods(image, r, c, nsize)
    neighborhood = [];
    n = (nsize-1)/2;
    newr=[];
    newc=[];
    
    for i = 1:length(r)
        if( r(i)>n && r(i)<size(image,1)-n && c(i)>n && c(i)<size(image,2)-n )
            nMatrix = image(r(i)-n:r(i)+n,c(i)-n:c(i)+n);
            nVector = reshape(nMatrix,1,[]);
            neighborhood = [ neighborhood; nVector];
            newr=[newr; r(i)];
            newc=[newc; c(i)];
        end
    end
    
    