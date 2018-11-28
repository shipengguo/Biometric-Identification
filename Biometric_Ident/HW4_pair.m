function [list,parameters] = minutiae_pairing(data1,data2,t_d,t_theta)
a = 1;

theta_T = data1(:,3) * 11.25;
theta_Q = data2(:,3) * 11.25;
x_T = data1(:,1);
x_Q = data2(:,1);
y_T = data1(:,2);
y_Q = data2(:,2);
stored = zeros(length(data1)*length(data2),5);

for i = 1:1:length(data1)
    for j = 1:1:length(data2)
        delta_theta = theta_T(i) - theta_Q(j);
        if delta_theta < 0
            delta_theta = delta_theta + 360;
        end
        delta_x = x_T(i) - x_Q(j)*cosd(delta_theta) - y_Q(j)*sind(delta_theta);
        delta_y = y_T(i) + x_Q(j)*sind(delta_theta) - y_Q(j)*cosd(delta_theta);
        stored(a,:) = [i,j,delta_x,delta_y,delta_theta];
        a = a + 1;
    end
end

x_interval = linspace(-300,300,41);
y_interval = linspace(-300,300,41);
theta_interval = linspace(0,360,21);
for i = 1:1:length(x_interval)-1
    for j = 1:1:length(y_interval)-1
        for k = 1:1:length(theta_interval)-1
            A(i,j,k) = 0;
        end
    end
end

history = cell(length(x_interval)-1,length(y_interval)-1,length(theta_interval)-1);
for i = 1:1:length(stored)
    for l = 2:1:length(theta_interval)
        if stored(i,5) >= theta_interval(l - 1) && stored(i,5) < theta_interval(l)
            for m = 2:1:length(x_interval)
                if stored(i,3) >= x_interval(m - 1) && stored(i,3) < x_interval(m)
                    for n = 2:1:length(y_interval)
                        if stored(i,4) >= y_interval(n - 1) && stored(i,4) < y_interval(n)
                           A(m-1,n-1,l-1) = A(m-1,n-1,l-1) + 1;
                           history{m-1,n-1,l-1} = [history{m-1,n-1,l-1};stored(i,3:5)];
                        end
                    end
                end
            end
        end
    end
end

[~,indmax]=max(A(:));
[m_max,n_max,j_max]=ind2sub(size(A),indmax);
ind = [m_max,n_max,j_max];

parameters = mean(cell2mat(history(ind(1),ind(2),ind(3))));
dx = parameters(1);
dy = parameters(2);
dtheta = parameters(3);


f_T = zeros(length(data1),1);
f_Q = zeros(length(data2),1);
count = 1;
k = 1;
delta_x = parameters(1);
delta_y = parameters(2);
delta_theta = parameters(3);

data2_rotated(:,1) = delta_x + data2(:,1)*cosd(delta_theta) + data2(:,2)*sind(delta_theta);
data2_rotated(:,2) = delta_y - data2(:,1)*sind(delta_theta) + data2(:,2)*cosd(delta_theta);
data2_rotated(:,3) = data2(:,3) + delta_theta;

list = [0,0,0,0];
for i = 1:1:length(data1)
    for j = 1:1:length(data2_rotated)
        distance = norm(data1(i,1:2)-data2_rotated(j,1:2));
        rotation = data1(i,3) - data2_rotated(j,3);
        if f_T(i) == 0 && f_Q(j) == 0 && distance <= t_d && abs(rotation) <= t_theta
            f_T(i) = 1;
            f_Q(j) = 1;
            list(count,:) = [i,j,distance,rotation];
            count = count + 1;
        end
        
        k = k + 1;
    end
end


end
