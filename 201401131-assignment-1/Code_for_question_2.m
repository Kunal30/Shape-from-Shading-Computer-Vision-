ps_list=[0.8389;0.5773;0.3638;0.1763;-0.1763;-0.3638;-0.5773;-0.8389];
qs_list=[0.7193;0.6363;0.5865;0.5596;-0.5596;-0.5865;-0.6363;-0.7193];
for i = 1:8
    ps_list(i)=ps_list(i)/sqrt(ps_list(i)^2+qs_list(i)^2+1);
    qs_list(i)=qs_list(i)/sqrt(ps_list(i)^2+qs_list(i)^2+1);
end
r = 24;
gridsize = 64;
mat = zeros(gridsize, gridsize);
mid = gridsize/2 ;

E=zeros(gridsize,gridsize,8);

for z=1:8
    
    for x = -31:32
        for y = -31:32
            if ((x)^2 + (y)^2) <= r^2 
                numerator = (((x)*ps_list(z)) + ((y) * qs_list(z)) + sqrt(r^2 - (x)^2 - (y)^2));
                denominator =  1.0 * (r * (1 + ps_list(z)^2 + qs_list(z)^2));
            
                answer = numerator/denominator;
                if answer<=0
                    answer = 0;
                end
                mat((mid+x),(mid+y)) = answer;
            end
        end
    end
   figure(z), imshow(mat)
    
    E(:,:,z)=mat;
end    

estimated_p=zeros(gridsize,gridsize);
estimated_q=zeros(gridsize,gridsize);
estimated_albedo=zeros(gridsize,gridsize);


source_mat=zeros(8,3);
source_mat=[-1.*(ps_list),-1.*(qs_list),ones(8,1)];


for i = 1:8
    source_mat(i,3)=source_mat(i,3)/sqrt(source_mat(i,1)^2+source_mat(i,2)^2+1);
end

disp(source_mat)
   
    for x=-31:32
        for y= -31:32
                    temp1 = reshape(E(mid+x,mid+y,:), [8,1]);
                    temp=pinv(source_mat)*(temp1);
                    
           
                    estimated_p(mid+x,mid+y)=-1*temp(1)/temp(3);
                    estimated_q(mid+x,mid+y)=-1*temp(2)/temp(3);
                    estimated_albedo(mid+x,mid+y)=temp(3)/sqrt(estimated_p((mid+x),(mid+y))^2+estimated_q((mid+x),(mid+y))^2+1);
                     if isnan(estimated_p(mid+x,mid+y))
                        estimated_p(mid+x,mid+y)=0;
                    end
                    if isnan(estimated_q(mid+x,mid+y))
                        estimated_q(mid+x,mid+y)=0;
                    end
                    if isnan(estimated_albedo(mid+x,mid+y))
                        estimated_albedo(mid+x,mid+y)=0;
                    end
        end
    end 
 
                
       
  disp(estimated_p)
  figure(9), imshow(estimated_p)
  figure(10), imshow(estimated_q)
  figure(11), imshow(estimated_albedo)

   


% SVD

    for x=-31:32
        for y= -31:32
                    temp1 = reshape(E(mid+x,mid+y,:), [8,1]);
                    [U S V] = svd(source_mat, 0);                                                           
                    temp= V *inv(S)*(U'*(temp1));                                        
                    estimated_p(mid+x,mid+y)=-1*temp(1)/temp(3);
                    estimated_q(mid+x,mid+y)=-1*temp(2)/temp(3);
                    estimated_albedo(mid+x,mid+y)=temp(3)/sqrt(estimated_p((mid+x),(mid+y))^2+estimated_q((mid+x),(mid+y))^2+1);
                    if isnan(estimated_p(mid+x,mid+y))
                        estimated_p(mid+x,mid+y)=0;
                    end
                    if isnan(estimated_q(mid+x,mid+y))
                        estimated_q(mid+x,mid+y)=0;
                    end
                    if isnan(estimated_albedo(mid+x,mid+y))
                        estimated_albedo(mid+x,mid+y)=0;
                    end
              
        end
    end 
    
    figure(12), imshow(estimated_p)
    figure(13), imshow(estimated_q)
    figure(14), imshow(estimated_albedo)


     
    % errors
    true_q = zeros(gridsize, gridsize);
    true_p = zeros(gridsize, gridsize);
    for x=-31:32
        for y= -31:32
            if(x^2+y^2<=r^2)
                z = sqrt(r^2 - x^2 - y^2);
            end
            if sqrt(r^2 - (x)^2 - (y)^2) >0
                true_p(mid + x, mid + y) = ( x / sqrt(r^2 - (x)^2 - (y)^2)) ;
                true_q(mid + x, mid + y) = ( y / sqrt(r^2 - (x)^2 - (y)^2)) ;
            end
        end
    end
 disp(true_p)
 figure(15),imshow(true_p);
 figure(16),imshow(true_q);
   
 
   P_Error = abs(true_p-estimated_p).^2;
   MSE_p = sum(P_Error(:))/numel(true_p);
   disp(MSE_p)
   
   Q_Error = abs(true_q-estimated_q).^2;
   MSE_q = sum(Q_Error(:))/numel(true_q);
   disp(MSE_q)