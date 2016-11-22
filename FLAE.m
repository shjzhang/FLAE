function [Q,W] = FLAE( r_base, b_base, weights )
    
    MM=zeros(3,3);
   
    for i=1:length(b_base(1,:))
        
        MM=MM+weights(i)*r_base(:,i)*b_base(:,i)';
    end
    
    Hx1=MM(1,1);    Hx2=MM(1,2);    Hx3=MM(1,3);
    Hy1=MM(2,1);    Hy2=MM(2,2);    Hy3=MM(2,3);
    Hz1=MM(3,1);    Hz2=MM(3,2);    Hz3=MM(3,3);
    
    W=[Hx1+Hy2+Hz3,-Hy3+Hz2,-Hz1+Hx3,-Hx2+Hy1;
       -Hy3+Hz2, Hx1-Hy2-Hz3,Hx2+Hy1,Hx3+Hz1;
       -Hz1+Hx3,Hx2+Hy1,Hy2-Hx1-Hz3,Hy3+Hz2;
       -Hx2+Hy1,Hx3+Hz1,Hy3+Hz2,Hz3-Hy2-Hx1];

   c=det(W);
   b=8*Hx3*Hy2*Hz1 - 8*Hx2*Hy3*Hz1 - 8*Hx3*Hy1*Hz2 + 8*Hx1*Hy3*Hz2 + 8*Hx2*Hy1*Hz3 - 8*Hx1*Hy2*Hz3;
   a=-2*Hx1*Hx1 - 2*Hx2*Hx2 - 2*Hx3*Hx3 - 2*Hy1*Hy1 - 2*Hy2*Hy2 - 2*Hy3*Hy3 - 2*Hz1*Hz1 - 2*Hz2*Hz2 - 2*Hz3*Hz3;

   temp1=2*a*a*a+27*b*b-72*a*c;
   DD=(real(sqrt(-4*(a*a+12*c))^3+temp1*temp1)+temp1)^(0.333333333333333);
   temp2=2.51984209978975*(a*a+12*c)/DD+1.5874010519682*DD;
   temp3=sqrt(-4*a+temp2);
   lambda=real(0.204124145231932*(temp3+real(sqrt(-8*a-temp2-29.3938769133981*b/temp3))));
   
   G=W-lambda*eye(4);
   
   
   
   pivot = G(1,1);  
   G(1,:) = G(1,:)/pivot;
   G(2,:) = G(2,:) - G(2, 1)*G(1,:);
   G(3,:) = G(3,:) - G(3, 1)*G(1,:);
   G(4,:) = G(4,:) - G(4, 1)*G(1,:);

   
   pivot = G(2,2);
   G(2,:) = G(2,:)/pivot;
   G(1,:) = G(1,:) - G(1, 2)*G(2,:);
   G(3,:) = G(3,:) - G(3, 2)*G(2,:);
   G(4,:) = G(4,:) - G(4, 2)*G(2,:);
   
   pivot = G(3,3);
   G(3,:) = G(3,:)/pivot;
   G(1,:) = G(1,:) - G(1, 3)*G(3,:);
   G(2,:) = G(2,:) - G(2, 3)*G(3,:);
   G(4,:) = G(4,:) - G(4, 3)*G(3,:);
   
   q=[G(1,4);G(2,4);G(3,4);-1];
 
   Q=q./norm(q);
   
   
   
end

