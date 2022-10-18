% This code consists of CCR profile comparison of experimental vs
% analytical having a converging cable routing

%% Defining Variables
E = 3500*10^6;
I = (pi*(3e-3)^4)/32;
L = 0.18;
z = 0:0.01:L;
F = 1.0;
b = 0.01;
a6 = -b/L;
r = 0.01;
% tau = 4;
tau = 40;

%% Governing equations for convergent profile
funx = @(s) cos((1/2).*(-b/L).*tau.*s.*s + L.*tau.*s);
funy = @(s) -sin((1/2).*(-b/L).*tau.*s.*s + L.*tau.*s);
x_value = zeros(1, length(z));
y_value = zeros(1, length(z));

for i = 1:length(z)
    x_value(i) = integral(funx, 0, z(i));
    y_value(i) = integral(funy, 0, z(i));
end

%% Image Processing
img = imread('C:\Users\C.ASHWIN\Desktop\flexible_robotics\Robot_images_testing\converging_profile\converging_3.jpg');
rot_img = imrotate(img, 360);
imshow(rot_img)
axis on
scaling_factor = 6000;
x_value_scaled = x_value*scaling_factor;
y_value_scaled = y_value*scaling_factor;
[input_x1, input_y1] = ginput(1);
hold on
plot(x_value_scaled-x_value_scaled(1)+input_x1, -y_value_scaled+y_value_scaled(1)+input_y1,'red')
hold off