clear,clc
close all

%load gallery and probe image sets
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
for i = 1 : length(srcFiles_2)
    filename = strcat('/Users/shipengguo/Desktop/Biometric/HW3/probeset/',srcFiles_2(i).name);
    probeset{i} =imread(filename);
end
probe=zeros(2500,200);
for i = 1:1:length(probeset)
    matrix = probeset{i};
    probe(:,i) = matrix(:);
end
score=zeros(200,100);
for i = 1:size(probe,2)
    for j = 1:size(gallery,2)
        score(i,j) = corr2(probe(:,i),gallery(:,j));
    end
end
% doing score level fusion
score_fusion=zeros(100,100);
for i=1:100
    for j=1:100
        if score(2*i-1,j)>score(2*i,j)
            score_fusion(i,j)=score(2*i-1,j);
        else
            score_fusion(i,j)=score(2*i,j);
        end
    end
end

genuine=diag(score_fusion);
for i = 1:size(score_fusion,1)
    imposter(i,:) = [score_fusion(i,1:i-1),score_fusion(i,i+1:end)];
end
imposter = imposter(:);

%plot score distribution

figure(1)
[N_g,X_g] = hist(genuine,10);
N_g = N_g/length(genuine);
plot(X_g,N_g,'Linewidth',2)
hold on
[N_i,X_i] = hist(imposter,10);
N_i = N_i/length(imposter);
plot(X_i,N_i,'r','Linewidth',2)
xlabel('match score')
ylabel('Probability')
title('Genuine & Imposter Score Distribution')
legend('Genuine','Imposter')

%plot CMC curve
threshold = 1:100;
prob = zeros(size(threshold));
for i = 1:1:length(threshold)
    count = 0;
    for j = 1:1:size(score_fusion,1)
        [~,index] = sort(score_fusion(j,:),'descend');
        indices = index(1:i);
        for k = 1:i
            if j == indices(k)
                count = count + 1;
            end
        end
    end
    prob(i) = count / size(score_fusion,1) * 100;
end
figure(2)
plot(threshold,prob,'Linewidth',2)
xlabel('Rank (t)')
ylabel('Rank-t Identification Rate (%)')
title('CMC Curve(multi-sample)')
axis([1 100 80 100]);
grid on;

%plot ROC curve
threshold = 0.2:0.001:1;
FRR = zeros(size(threshold));
for j = 1:1:length(threshold)
    count = 0;
    for k = 1:1:length(genuine)
        if genuine(k) < threshold(j)
            count = count + 1;
        end
    end
    FRR(j) = count / length(genuine);
end

FAR = zeros(size(threshold));
for j = 1:1:length(threshold)
    count = 0;
    for k = 1:1:length(imposter)
        if imposter(k) >= threshold(j)
            count = count + 1;
        end
    end
    FAR(j) = count / length(imposter);
end
figure(3)
plot(FAR,FRR,'Linewidth',2)
hold on
plot(0:0.001:0.5,0:.001:0.5,'-.r')
legend('ROC Curve','y = x')
xlabel('False Accept Rate (FAR)')
ylabel('False Reject Rate (FRR)')
title('ROC Curve(multi-sample)')
grid on;