clear,clc
close all

%face recognition
srcFiles_1 = dir('/Users/shipengguo/Desktop/Biometric/HW3/galleryset/*.pgm');
galleryset=cell(length(srcFiles_1),1);
for i = 1 : length(srcFiles_1)
    filename = strcat('/Users/shipengguo/Desktop/Biometric/HW3/galleryset/',srcFiles_1(i).name);
    galleryset{i} = imread(filename);
end
gallery=zeros(2500,100);
for i = 1:1:length(galleryset)
    matrix = galleryset{i};
    gallery(:,i) =matrix(:);
end

srcFiles_2 = dir('/Users/shipengguo/Desktop/Biometric/HW3/probeset/*.pgm');
probeset=cell(length(srcFiles_2),1);
for i = 1:length(srcFiles_2)
    filename = strcat('/Users/shipengguo/Desktop/Biometric/HW3/probeset/',srcFiles_2(i).name);
    probeset{i} =imread(filename);
end
probe=zeros(2500,100);
for i = 1:100
    matrix = probeset{2*i-1};
    probe(:,i) = matrix(:);
end
%caculate the corr score
score_corr=zeros(100,100);
for i = 1:100
    for j = 1:size(gallery,2)
        score_corr(i,j) = corr2(probe(:,i),gallery(:,j));
    end
end

%iris recognition
clc;
clear all;

srcFiles_3 = dir('/Users/shipengguo/Desktop/Biometric/HW5/iris_img/gallery/*.tiff');
gallery_enroll=cell(length(srcFiles_3),2);
for i = 1 : length(srcFiles_3)
    eyeimage_gallery = strcat('/Users/shipengguo/Desktop/Biometric/HW5/iris_img/gallery/',srcFiles_3(i).name);
    [template,mask]=createiristemplate(eyeimage_gallery);
    gallery_enroll{i,1}=template;
    gallery_enroll{i,2}=mask;
    
end

srcFiles_4=dir('/Users/shipengguo/Desktop/Biometric/HW5/iris_img/probe/*.tiff');
probe_enroll=cell(length(srcFiles_4),2);
for i=1:length(srcFiles_4)
    eyeimage_probe = strcat('/Users/shipengguo/Desktop/Biometric/HW5/iris_img/probe/',srcFiles_4(i).name);
    [template,mask]=createiristemplate(eyeimage_probe);
    probe_enroll{i,1}=template;
    probe_enroll{i,2}=mask;
end

hd=zeros(length(srcFiles_4),length(srcFiles_3));
for i=1:length(srcFiles_4)
    for j=1:length(srcFiles_3)
        hd(i,j)=gethammingdistance(probe_enroll{i,1},probe_enroll{i,2},gallery_enroll{j,1},gallery_enroll{j,2},1);
        
    end
end
genuine=zeros(length(hd),1);
genuine(:) = diag(hd);
for i = 1:1:size(hd,1)
    imposter(i,:) = [hd(i,1:i-1),hd(i,i+1:end)];
end
imposter = imposter(:);
mean_gen=mean(genuine);
mean_imp=mean(imposter);
std_gen=std(genuine);
std_imp=std(imposter);
figure(2)
x_1=linspace(min(genuine),max(genuine));
y_1=normpdf(x_1,mean_gen,std_gen);%/length(genuine);
plot(x_1,y_1,'LineWidth',2)
hold on;
x_2=linspace(min(imposter),max(imposter));
y_2=normpdf(x_2,mean_imp,std_imp);%/length(imposter);
plot(x_2,y_2,'r','LineWidth',2)
xlabel('hamming distance')
ylabel('probability density')
legend('genuine','imposter')