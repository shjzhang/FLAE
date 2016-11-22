function [q,K,iter,time_poly,time_eig]=QUEST_shuster(sxy,cxy,weight)

    tic;
    B=0;
    nn=length(cxy(1,:));
    for j=1:nn
        B=B+weight(j)*sxy(:,j)*cxy(:,j)';
    end

    S=B'+B;
    sigma=B(1,1)+B(2,2)+B(3,3);
    Z=[B(2,3)-B(3,2),B(3,1)-B(1,3),B(1,2)-B(2,1)]';
    adjS=det(S)*inv(S);
    ka=trace(adjS);
    delta=det(S);
    K=[S-sigma*eye(3,3),Z;Z',sigma];

    a=sigma*sigma-ka;
    b=sigma*sigma+Z'*Z;
    c=delta+Z'*S*Z;
    d=Z'*S*S*Z;

    time_poly=toc;
    
    tic;
    lambda=1.0;
    old_lambda=0.0;
    iter=0;

    while(abs(old_lambda-lambda)>1e-5)
        old_lambda=lambda;
        lambda=lambda-((lambda^4-(a+b)*lambda^2-c*lambda+(a*b+c*sigma-d))/(4*lambda^3-2*(a+b)*lambda-c));
        iter=iter+1;
    end


    arf=lambda^2-sigma^2+ka;
    bet=lambda-sigma;
    X=(arf*eye(3,3)+bet*S+S^2)*Z;
    miu=(lambda+sigma)*arf-det(S);
    q=[miu;X]/sqrt(miu^2+X(1)^2+X(2)^2+X(3)^2);
    time_eig=toc;

end