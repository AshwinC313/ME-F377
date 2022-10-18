%% 0) Parallel Tendon deflections
E = 3500*10^6;
I = (pi*(3e-3)^4)/32;
% I1 = (pi*(3e-3)^4)/32;
L = 0.6;
z = 0:0.01:L;
F = 1.0;
b = 0.035;
a6 = -b/L;
r = 0.01;
% tau = 4;
tau = 5.85;
% delta = (F*z.^3)/(3*E*I);
% delta1 = delta - (z/L).*((F*z.^3)/(4*E*I));
% dc = (12*E*I/(L^3))*(delta/F);
% dc1 = (12*E*I/(L^3))*(delta1/F);
% dll = z/L;
% figure
% plot(dll, dc)
% hold on
% plot(dll, dc1)
% xlim([0 1]);
% ylim([0 1]);
% hold off
%% 1) Non Parallel Tendon deflections
% a1 = 0;
% a2 = -0.2*(b/L);
% a3 = -0.4*(b/L);
% a4 = -0.6*(b/L);
% a5 = -0.8*(b/L);


% delta_non_p1 = calc_delta_non_parallel(a1);
% delta_non_p2 = calc_delta_non_parallel(a2);
% delta_non_p3 = calc_delta_non_parallel(a3);
% delta_non_p4 = calc_delta_non_parallel(a4);
% delta_non_p5 = calc_delta_non_parallel(a5);
% delta_non_p6 = calc_delta_non_parallel(a6);
% 
% dc_np1 = (12*E*I/(L^3))*(delta_non_p1/F);
% dc_np2 = (12*E*I/(L^3))*(delta_non_p2/F);
% dc_np3 = (12*E*I/(L^3))*(delta_non_p3/F);
% dc_np4 = (12*E*I/(L^3))*(delta_non_p4/F);
% dc_np5 = (12*E*I/(L^3))*(delta_non_p5/F);
% dc_np6 = (12*E*I/(L^3))*(delta_non_p6/F);

% figure
% plot(dll, dc_np1);
% hold on
% plot(dll, dc_np2);
% plot(dll, dc_np3);
% plot(dll, dc_np4);
% plot(dll, dc_np5);
% plot(dll, dc_np6);
% xlim([0 1]);
% ylim([0 1]);
% legend('a=0 (parallel)', 'a = -0.2b/L', 'a = -0.4b/L', 'a = -0.6b/L', 'a=-0.8b/L','a=-b/L(converging)'); 
% hold off

%% 2) relationship between a tendon’s stiffness and the deflection response to a tip load for both parallel routing and converging routing
% r = 0.1:0.01:1;
% b = 0.1:0.01:1;
% k = 2.3*10e6;
% delta_parallel_tip = (F*L^3)/(3*E*I) - ((r.*r.*k*L)./(E*I + r.*r.*k*L))*((F*L^3)/(4*E*I));
% delta_converging_tip = (F*L^3)/(3*E*I) - ((b.*b.*k*L)./(3*E*I + b.*b.*k*L))*((F*L^3)/(3*E*I));
% dtc_parallel = ((12*E*I)/(L^3))*(delta_parallel_tip/F);
% dtc_converging = ((12*E*I)/(L^3))*(delta_converging_tip/F);
% dts = ((r.^2)*k*L)/(E*I);
% figure
% plot(dts, dtc_parallel);
% hold on
% plot(dts, dtc_converging);
% xlim([0 50]);
% ylim([0 4]);
% legend('parallel','converging');
% hold off

%% 3) Profile of flexible rod for parallel tendons
% solving quadratic equation of theta in order to get relation of theta
% along x and y axis

% r = 0.01;
% tau = 10;
% % quadratic equation of theta
% alpha = (F*r*r)/(2*E*I);
% beta = 1 - (pi*F*r*r)/(2*E*I) + r*r*tau - r*F*z;
% gamma = F*r*r*pi*pi/(8*E*I) - (pi*r*r*tau/2) + (r*F*z/2);
% theta_sol = (-beta + sqrt(beta.*beta - 4*alpha*gamma))./(2*alpha);
% figure
% hold on
% plot(r*cos(theta_sol), r*sin(theta_sol));
% title('profile of the backbone');
% hold off
%% 4) Modified profile of flexible rod for parallel tendons
x_s = ((E*I)/(r*tau))*sin((r*tau*z)/(E*I));
y_s = ((E*I)/(r*tau))*cos((r*tau*z)/(E*I));
coordinate_matrix_parallel_routing = [x_s ;y_s];
save example.mat coordinate_matrix_parallel_routing -v7.3;
% x_s1 = ((E*I)/(r*tau))*cos((r*tau*z)/(E*I));
% y_s1 = ((E*I)/(r*tau))*sin((r*tau*z)/(E*I));
% figure
% hold on
% plot(x_s, y_s);
% title('modified profile of the robot');
% hold off

%% 5) Profile of flexible rod for converging tendons
B = L;
funx = @(s) cos(0.5.*a6.*tau.*s.*s + B.*tau.*s);
funy = @(s) sin(0.5.*a6.*tau.*s.*s + B.*tau.*s);
integralx_storage = zeros(1, length(z));
integraly_storage = zeros(1, length(z));
x_s2 = zeros(1, length(z));
y_s2 = zeros(1, length(z));
for i = 1:length(z)
    integralx_storage(i) = integral(funx, 0, z(i));
    integraly_storage(i) = integral(funy, 0, z(i));
    x_s2(i) = integral(funy, 0, z(i));
    y_s2(i) = integral(funx, 0, z(i));
end

% figure
% hold on
% title('Profile of the rod having converging tendons');
% plot(integralx_storage, integraly_storage);
% hold off

% subplot(2,1,1)
% plot(x_s, y_s);
% % plot(x_s1, y_s1);
% plot(y_s, x_s);
% title('Profile of the robot - parallel tendon routing');
% axis equal
% subplot(2,1,2);
% % plot(integralx_storage, integraly_storage);
% plot(x_s2, y_s2);
% % plot(integraly_storage, integralx_storage);
% title('Profile of the robot - converging tendon routing');
% hold off

%% 6) Experimental comparison of analytically determined curves with the actual robot profile
img1 = imread('C:\Users\C.ASHWIN\Desktop\flexible_robotics\Robot_images_testing\converging_profile\converging_3.jpg');
J1 = imrotate(img1,360);
% imshow(img)
imshow(J1)
% J2 = imresize(J1, 0.05);
% imshow(J2)
axis on
scaling_factor = 2000;
x_s_mod = x_s*scaling_factor;
y_s_mod = y_s*scaling_factor;
[x,y] = ginput(1);

% for i = 1:length(x_s_mod)
%     x_s_mod(i) = x_s_mod(i) + x;
%     y_s_mod(i) = y_s_mod(i) - y_s_mod(1) + y;
% end
hold on
plot(x_s_mod-x_s_mod(1)+x, -y_s_mod+y_s_mod(1)+y,'red');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [X,Y] = ginput(3);
% [a,b,c] = polyn_predictor(X,Y);
% y_exp = a*z.*z + b*z + c;
% 
% plot(x_s, y_exp - c);
% % plot(x_s2, y_exp - c);
% hold on
% plot(x_s, y_s - y_s(1));
% % plot(x_s2, y_s2 - y_s2(1));
% hold off
% legend('robot profile obtained experimentally', 'robot profile obtained analytically');
% xlabel('x coordinate of robot(in m)');
% ylabel('y coordinate of robot(in m)');

%% 7) Necessary functions used
function delta_non_parallel = calc_delta_non_parallel(a)
    E = 3500*10^6;
    I = 10^-5;
    L = 0.095;
    z = 0:0.001:0.15;
    F = 1.0;
    b = 0.01;
    delta_non_parallel = (F*z.^3)/(3*E*I) - (((a*z + 3*b).^2).*z)./(3*L*(a*a*L*L + 3*a*b*b + 3*b*b)).*(F*z.^3)/(4*E*I);
end

% This function returns the coefficients of the predicted polynomial
% assumption before creating this function : robot profile observed in the
% figure is assumed to be following a quadratic trajectory.
function [a, b, c] = polyn_predictor(x,y)
    % input : x and y arrays containing values of the coordinates of the
    % polynomial (x1, x2, x3) and (y1, y2, y3)
    a = (1/(x(1) - x(3)))*(((y(1) - y(2))/(x(1) - x(2))) - ((y(2) - y(3))/(x(2) - x(3))));
    b = ((y(1) - y(2))/(x(1) - x(2))) - a*(x(1) + x(2));
    c = y(1) - a*x(1)*x(1) - b*x(1);
end
