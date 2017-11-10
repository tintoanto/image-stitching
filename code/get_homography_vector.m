function h = get_homography_vector(pm, pm_random4idx)

    im1x = pm{1}(:,1);
    im1y = pm{1}(:,2);
    x1 = im1x(pm_random4idx);
    y1 = im1y(pm_random4idx);

    im2x = pm{2}(:,1);
    im2y = pm{2}(:,2);
    x2 = im2x(pm_random4idx);
    y2 = im2y(pm_random4idx);

    A = [];
    B = [];
    for i = 1:4     
        A_2 = [ x1(i)   y1(i)   1   0       0       0   -(x1(i)*x2(i))  -(y1(i)*x2(i)); ...
                0       0       0   x1(i)   y1(i)   1   -(x1(i)*y2(i))  -(y1(i)*y2(i))    ];
        A   = [ A ; A_2];
        
        B_2 = [ x2(i)   ;   y2(i)   ];
        B   = [ B       ;   B_2     ];
        
    end
    
    hv = A\B;
    h = [   hv(1)   hv(2)   hv(3)   ;   
            hv(4)   hv(5)   hv(6)   ;   
            hv(7)   hv(8)   1       ]; 
end