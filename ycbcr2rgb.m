function RGB = ycbcr2rgb(Cb,Cr)
    if ~isequal(size(Cb),size(Cr))
        error('Cb and Cr must have the same size.')
    end
    len = length(Cb(:));
    
    if min(Cb,[],'all') < 0
        error('Cb must be >= 0.')
    end
    if max(Cb,[],'all') > 1
        error('Cb must be <= 1.')
    end
    Cb = Cb*256;
    if min(Cr,[],'all') < 0
        error('Cr must be >= 0.')
    end
    if max(Cr,[],'all') > 1
        error('Cr must be <= 1.')
    end
    Cr = Cr*256;
    
    Y = 256;
    A = [   65.783,     129.057,    25.064;...
            -37.945,    -74.494,    112.439;...
            112.439,    -94.154,    -18.285];
    R = NaN(size(Cb));
    G = NaN(size(Cb));
    B = NaN(size(Cb));
    
    for i=1:len
        ycbcr = [Y;Cb(i);Cr(i)] - [16;128;128];
        rgb = linsolve(A/256,ycbcr);
        
        R(i) = rgb(1);
        G(i) = rgb(2);
        B(i) = rgb(3);
    end
    
    RGB(:,:,1) = R/512;
    RGB(:,:,2) = G/512;
    RGB(:,:,3) = B/512;

end