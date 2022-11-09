% This code consists of CCR profile comparison of experimental vs
% analytical having a converging cable routing

%% Defining Variables
E = 3500*10^6;
I = (pi*(3e-3)^4)/32;
L = 0.18;
z = 0:0.01:L;
b = 0.005;
a = b*2;
tau = 55;
%% Governing equations for convergent profile
funx = @(s) cos((-(a-b).*tau.*s.*s)/(2.*L.*E.*I) - (b.*tau.*s)/(E.*I));
funy = @(s) -sin((-(a-b).*tau.*s.*s)/(2.*L.*E.*I) - (b.*tau.*s)/(E.*I));
x_value = zeros(1, length(z));
y_value = zeros(1, length(z));

for i = 1:length(z)
    x_value(i) = integral(funx, 0, z(i));
    y_value(i) = integral(funy, 0, z(i));
end
%% Image Processing
img = imread('C:\Users\C.ASHWIN\Desktop\flexible_robotics\Robot_images_testing\diverging_profile\diverging_2.jpg');
% rot_img = imrotate(img, 360);
% imshow(rot_img)
imshow(img)
axis on
scaling_factor = 7000;
x_value_scaled = x_value*scaling_factor;
y_value_scaled = y_value*scaling_factor;
[input_x1, input_y1] = ginput(1);
hold on
plot(x_value_scaled-x_value_scaled(1)+input_x1, -y_value_scaled+y_value_scaled(1)+input_y1,'red')
hold off