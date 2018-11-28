clc;
clear all;

%part 2 Hough algorithm

x1=load('FP(1).txt');
x2=load('FP(6).txt');

%处理角度
x1(:,3)=x1(:,3)*11.25;
x2(:,3)=x2(:,3)*11.25;

delta_theta=[];
delta_x=[];
delta_y=[];

%得到所有的△x,△y和△theta组合
for i=1:length(x1)
    for j=1:length(x2)
        theta=x1(i,3)-x2(j,3);
        delta_theta=[delta_theta;theta];
        x=x1(i,1)-x2(j,1)*cos(theta*pi/180)-x2(j,2)*sin(theta*pi/180);
        delta_x=[delta_x;x];
        y=x1(i,2)+x2(j,1)*sin(theta*pi/180)-x2(j,2)*cos(theta*pi/180);
        delta_y=[delta_y;y];
    end
end

%△theta转到（0，360）
for i=1:length(delta_theta)
    if delta_theta(i)<0
        delta_theta(i)=delta_theta(i)+360;
    end
end
interval_x=10;
interval_y=10;
count=zeros(1000/interval_x,1000/interval_y,33);
%count_theta=zeros(17,1);
store=cell(1000/interval_x,1000/interval_y,33);

for i=1:length(x1)*length(x2)
    for j=1:33
         if delta_theta(i)>=(-11.25+11.25*j)&&delta_theta(i)<(11.25*j)
             %count_theta(j)=count_theta(j)+1;
             for m=1:1000/interval_x
                 if delta_x(i)>=(-400+m*interval_x)&&delta_x(i)<(-380+m*interval_x)
                     for n=1:1000/interval_y
                         if delta_y(i)>=(-400+n*interval_y)&&delta_y(i)<(-380+n*interval_y)
                            count(m,n,j)=count(m,n,j)+1;
                            store{m,n,j}=[store{m,n,j};[delta_x(i) delta_y(i) delta_theta(i)]];
                         end
                     end
                 end
             end
         end
    end
end

%最大的bin的计数数值
count_max=max(max(max(count)))
%找最大值的坐标
[count_max,indmax]=max(count(:));

[m_max,n_max,j_max]=ind2sub(size(count),indmax)
% 输出这个最大计数bin的所有sets
store{m_max,n_max,j_max}
% 根据这些sets取hough algorithm的最终输出△x,△y和△theta,取均值
final_val=[];
final_val=[final_val;mean(store{m_max,n_max,j_max},1)]



%part 3 paring
%先alignment
x2_new=zeros(length(x2),3);

x2_new(:,3)=x2(:,3)+final_val(1,3);
x2_new(:,1)=final_val(1,1)+x2(:,1)*cosd(final_val(1,3))+x2(:,2)*sind(final_val(1,3));
x2_new(:,2)=final_val(1,2)+x2(:,2)*cosd(final_val(1,3))-x2(:,1)*sind(final_val(1,3));
    
%find pairs
distance=zeros(length(x1),length(x2));
list=[];
f_T=zeros(length(x1));
f_Q=zeros(length(x2));
%store_dist=[];
%ful_thrh=cell(length(x1),1);
for i=1:length(x2)
    for j=1:length(x1)
        distance(i,j)=norm(x2_new(i,1:2)-x1(j,1:2));
        if f_T(j)==0&&f_Q(i)==0&&distance(i,j)<=10&&abs(x1(j,3)-x2_new(i,3))<=12
            f_T(j)=1;
            f_Q(i)=1;
            list=[list;[i,j]];
            
            
    %[min_distance,min_ind]=min(ful_thrh(i));
    %list_final=[list_final;list(min_ind)];
    %x2_new(list(min_ind,2),:)=zeros(1,3);
    %count=count+1;
        end
    end
end

paired_points=length(list)*2
unpaired_points=length(x1)+length(x2)-paired_points


  
        
    
