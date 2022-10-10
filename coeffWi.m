function out = coeffWi(ii,a)
    W8 = fliplr([32768, -131072,212992, - 180224,84480, - 21504,2688, - 128]);
    W7 = fliplr([0,8192, - 28672,39424, - 26880,9408, - 1568,98]);
    W6 = fliplr([0,0,2048, - 6144,6912, - 3584,840, - 72]);
    W5 = fliplr([0,0,0,512, - 1280,1120, - 400,50]);
    W4 = fliplr([0,0,0,0,128, - 256,160, - 32]);
    W3 = fliplr([0,0,0,0,0,32, - 48,18]);
    W2 = fliplr([0,0,0,0,0,0,8, - 8]);
    W1 = fliplr([0,0,0,0,0,0,0,2]);
    W = [W1;W2;W3;W4;W5;W6;W7;W8];
    
    if (ii >= 1)
        tmp = 0;
        for q = ii:8
            tmp = tmp + besseli(q,a/2)*W(q,ii);
        end
        tmp = 2 * tmp;
    else
        tmp = 0;
        for q = 0:8
            if (q >= 1)
                tmp = tmp + 2*(-1)^q*besseli(q,a/2);
            else
                tmp = tmp + besseli(q,a/2);
            end
        end
    end
    out = tmp;
end