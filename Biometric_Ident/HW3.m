clear,clc
close all
%I_gal=cell(1,100);
%for i=1:100
%    img_path=strcat(['/Users/shipengguo/Desktop/Biometric/HW3/galleryset/subject',num2str(i),'_img1.pgm']);
%    I_gal{i}=imread(img_path);
%end
srcFiles_1 = dir('/Users/shipengguo/Desktop/Biometric/HW3/galleryset/*.pgm');
for i = 1 : length(srcFiles_1)
    filename = strcat('/Users/shipengguo/Desktop/Biometric/HW3/galleryset/',srcFiles_1(i).name);
    galleryset{i} = imread(filename);
end 
for i = 1:length(galleryset)
    matrix = galleryset{i};
    ma_top=matrix(:,:);     %更改输入的图片大小
    gallery(:,i) = ma_top(:);
end
gallery = double(gallery);
x_=mean(gallery,2);
for i=1:length(galleryset)
    
    gallery_(:,i) = gallery(:,i)-x_;%减去均值
end


[V D latent]=pca(gallery_');%转置是因为这里的函数用一个行向量代表一个图像
sum1=sum(latent);
sum_=0;
i_gallery=1;
while (sum_/sum1<0.99)
    sum_=sum_+latent(i_gallery,1);
    i_gallery=i_gallery+1;
end
i_gallery=i_gallery-1;
u = V(:,1:i_gallery);
for i=1:size(u,2)
    u(:,i)=u(:,i)/norm(u(:,i));
end
a_gallery = u'*gallery_;
%处理probe
srcFiles_2 = dir('/Users/shipengguo/Desktop/Biometric/HW3/probeset/*.pgm');
for i = 1 : length(srcFiles_2)
    filename = strcat('/Users/shipengguo/Desktop/Biometric/HW3/probeset/',srcFiles_2(i).name);
    probeset{i} = imread(filename);
end

for i = 1:length(probeset)
    matrix = probeset{i};
    ma_top=matrix(:,:); %更改图片大小
    probe(:,i) = ma_top(:);
end

x_probe=mean(probe,2);
probe = double(probe);
for i=1:length(probeset)
   
    probe_(:,i) = probe(:,i)-x_probe;
end

% [Vp Dp latent_p]=princomp(probe_');
% sum1=sum(latent_p);
% sum_=0;
% i_probe=1;
% while (sum_/sum1<0.9)
%     sum_=sum_+latent_p(i_probe,1);
%     i_probe=i_probe+1;
% end
% t=max(i_probe,i_gallery);
% u=[];
% u=[u V(:,1:t)];%gallery的将维度后选取的特征向量

%找对应的weight
% for i=1:length(I_gal)
%     a=u'*gallery_(:,i);
%     a_gallery(:,i)=a;
% end
% u_probe=u;

%找weight
% for i=1:length(probeset)
%     b=u_probe'*probe_(:,i);
%     b_prb(:,i)=b;
% end
b_prb = u'*probe_;
%用欧拉距离算出match score
for i=1:length(galleryset)
    for j=1:length(probeset)
        score(i,j)=norm(a_gallery(:,i)-b_prb(:,j));
    end
end



for i = 1:1:size(score,1)
    genuine(i,:) = score(i,2*i-1:2*i);
    imposter(i,:) = [score(i,1:2*i-2),score(i,2*i+1:end)];
end
genuine = genuine(:);
imposter = imposter(:);




% Match Score Distribution 
figure(1)
[N_g,X_g] = hist(genuine,10);
N_g = N_g/length(genuine);
N = 2;
plot(X_g,N_g,'Linewidth',1)
hold on
[N_i,X_i] = hist(imposter,10);
N_i = N_i/length(imposter);
plot(X_i,N_i,'Linewidth',1)
xlabel('Eu distance')
ylabel('Probability')
title('Score Distribution(right face)')
legend('Genuine','Imposter')

% CMC curve
score=score';

threshold = 1:100;
prob = zeros(size(threshold));
for i = 1:1:length(threshold)
    count = 0;
    for j = 1:1:size(score,1)
        [~,index] = sort(score(j,:),'ascend');
        indices = index(1:i);
        for k = 1:i
            if (j+1)/2 == indices(k)
                count = count + 1;
            end
            if j/2 == indices(k)
                count = count + 1;
            end
        end
    end
    prob(i) = count / size(score,1) * 100;
end
figure(2)
plot(threshold,prob,'Linewidth',1)
xlabel('Rank (t)')
ylabel('Rank-t Identification Rate (%)')
title('CMC Curve(right face)')
grid on

% ROC curve
threshold = 1000:50:4000;
FRR = zeros(size(threshold));
for j = 1:1:length(threshold)
    count = 0;
    for k = 1:1:length(genuine)
        if genuine(k) >= threshold(j)
            count = count + 1;
        end
    end
    FRR(j) = count / length(genuine);
    
end

FAR = zeros(size(threshold));
for j = 1:1:length(threshold)
    count = 0;
    for k = 1:1:length(imposter)
        if imposter(k) < threshold(j)
            count = count + 1;
        end
    end
    FAR(j) = count / length(imposter);
end
figure(3)
plot(FAR,FRR,'Linewidth',1)
hold on
plot(0:0.001:0.5,0:0.001:0.5,'-.r')
legend('ROC Curve','y = x')
xlabel('False Accept Rate (FAR)')
ylabel('False Reject Rate (FRR)')
title('ROC Curve(right face)')

%计算d'
u0=mean(genuine);
u1=mean(imposter);
d0=std(genuine);
d1=std(imposter);
d_prime=sqrt(2)*(u0-u1)/sqrt(d0^2+d1^2);

