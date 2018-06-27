r = 24;
gridsize = 64;
E = zeros(gridsize, gridsize);
mid = gridsize/2 ;
    
    for x = -31:32
        for y = -31:32
            if ((x)^2 + (y)^2) <= r^2 
                numerator = 1.0*( sqrt(r^2 - (x)^2 - (y)^2));
                denominator =  1.0*r;
            
                answer = numerator/denominator;
                if answer<=0
                    answer = 0;
                end
                E((mid+x),(mid+y)) = answer;
            end
        end
    end        
imshow(E);   
disp(E);
  