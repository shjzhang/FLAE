clear all
clc

addpath('quaternion_library'); 
ll=10000;


C_ture=[0.352, 0.864,0.360;
        -0.864,0.152,0.48;
        0.36,-0.48,0.8;];
    
r1=[1;0;0];
r2=[0;1;0];
r3=[0;0;1];

n1=wgn(ll,3,1,'linear')/220;
n2=wgn(ll,3,1,'linear')/220;
n3=wgn(ll,3,1,'linear')/220;


eig_time_QUEST=zeros(ll,1);
eig_time_FLAE=zeros(ll,1);

time_QUEST=zeros(ll,1);
time_FLAE=zeros(ll,1);

quaternion_QUEST=zeros(ll,4);
quaternion_FLAE=zeros(ll,4);

iter_QUEST=zeros(ll,1);
iter_FLAE=zeros(ll,1);

time_poly_QUEST=zeros(ll,1);
time_poly_FLAE=zeros(ll,1);

time_eig_QUEST=zeros(ll,1);
time_eig_FLAE=zeros(ll,1);

for i=1:ll
    
    b1=C_ture*r1+n1(i,:)';
    b2=C_ture*r2+n2(i,:)';
    b3=C_ture*r3+n3(i,:)';
    
    len=3;
    Db=[b1,b2,b3];
    Dr=[r1,r2,r3];
    
    a = 1.*rand([len,1])+0;
    a = a./sum(a);
    weights=a;
    
    tic;
    [q_QUEST,K,iter_QUEST(i),time_poly_QUEST(i),time_eig_QUEST(i)]=QUEST_shuster(Dr,Db,weights);
    time_QUEST(i)=toc;
    
    tic;
    [q_FLAE,W,iter_FLAE(i),time_poly_FLAE(i),time_eig_FLAE(i)]=FLAE_newton(Dr,Db,weights);
    time_FLAE(i)=toc;
    
    tic;
    eig(K);
    eig_time_QUEST(i)=toc;
    
    tic;
    eig(W);
    eig_time_FLAE(i)=toc;
    
    quaternion_QUEST(i,:)=q_QUEST';
    quaternion_FLAE(i,:)=q_FLAE';
    
end

euler_QUEST=quatern2euler(quaternConj(quaternion_QUEST))*180/pi;
euler_FLAE=quatern2euler(quaternConj(quaternion_FLAE))*180/pi;

figure(1);
plot(time_QUEST(100:ll),'*'); hold on;
plot(time_FLAE(100:ll),'o'); hold off;
legend('QUEST','FLAE');
ylim([0 1e-4]);
title('Overall Time');
ylabel('Time (s)');
xlabel('Sample Index');

figure(2);
plot(eig_time_QUEST(100:ll),'*'); hold on;
plot(eig_time_FLAE(100:ll),'o'); hold off;
legend('QUEST','FLAE');
ylim([0 1e-5]);
title('Time Consumption of eig(K) and eig(W)');
ylabel('Time (s)');
xlabel('Sample Index');

figure(3);
subplot(2,2,1);
plot(quaternion_QUEST);
title('Quaternion QUEST');
subplot(2,2,2);
plot(quaternion_FLAE);
title('Quaternion FLAE');
subplot(2,2,3);
plot(euler_QUEST);
title('Euler Angles QUEST');
subplot(2,2,4);
plot(euler_FLAE);
title('Euler Angles FLAE');

figure(4);
plot(time_poly_QUEST(100:ll),'*'); hold on;
plot(time_poly_FLAE(100:ll),'o'); hold off;
legend('QUEST','FLAE');
ylim([0 1e-4]);
title('Time Consumption for Computing Polynomial Coefficients');
ylabel('Time (s)');
xlabel('Sample Index');

figure(5);
plot(time_eig_QUEST(100:ll),'*'); hold on;
plot(time_eig_FLAE(100:ll),'o'); hold off;
legend('QUEST','FLAE');
ylim([0 1e-4]);
title('Time Consumption for Computing Newton Iter.');
ylabel('Time (s)');
xlabel('Sample Index');

figure(6);
subplot(2,1,1);
plot(iter_QUEST,'*'); hold off;
title('Iteration Amounts of QUEST Newton Algorithm (eps=1e-5)');
ylabel('Iteration Amount');
xlabel('Sample Index');

subplot(2,1,2);
plot(iter_FLAE,'o'); hold off;
title('Iteration Amounts of FLAE Newton Algorithm (eps=1e-5)');
ylabel('Iteration Amount');
xlabel('Sample Index');




mean(time_QUEST)
mean(time_FLAE)

mean(eig_time_QUEST)
mean(eig_time_FLAE)

    