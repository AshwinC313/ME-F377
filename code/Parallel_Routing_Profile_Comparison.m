% Profile comparison between the experimental vs analytical for the CCR
% having parallel routing.

%% Defining necessary variables
E = 3500*10^6; % Modulus of Elasticity
I = (pi*(3e-3)^4)/32; % Moment of Inertia
L = 0.6; % Length of CCR
z = 0:0.01:L; 
F = 1.0; % External load (if necessary) applied on the robot
b = 0.035; 
a6 = -b/L;
r = 0.01;
% tau = 4; % tension applied on the tendons 
tau = 5.85; % tension applied on the tendons
%% Relation of x with arc length parameter of CCR
x_s = ((E*I)/(r*tau))*sin((r*tau*z)/(E*I));
y_s = ((E*I)/(r*tau))*cos((r*tau*z)/(E*I));

%% Image processing for comparison
img = imread('C:\Users\C.ASHWIN\Desktop\flexible_robotics\Robot_images_testing\converging_profile\converging_3.jpg');
J1 = imrotate(img,360);
imshow(J1)
axis on
scaling_factor = 2000;
x_s_scaled = x_s*scaling_factor;
y_s_scaled = y_s*scaling_factor;
[input_x, input_y] = ginput(1);
hold on
plot(x_s_scaled - x_s_scaled(1) + input_x, -y_s_scaled + y_s_scaled(1) + input_y, 'red');
hold off