%% Parallel Tendon deflections
E = 3500*10^6;
I = 10^-5;
L = 0.15;
z = 0:0.001:0.15;
F = 1.0;
b = 0.035;
delta = (F*z.^3)/(3*E*I);
delta1 = delta - (z/L).*((F*z.^3)/(4*E*I));
dc = (12*E*I/(L^3))*(delta/F);
dc1 = (12*E*I/(L^3))*(delta1/F);
dll = z/L;
% figure
% plot(dll, dc)
% hold on
% plot(dll, dc1)
% xlim([0 1]);
% ylim([0 1]);
% hold off
%% Non Parallel Tendon deflections
a1 = 0;
a2 = -0.2*(b/L);
a3 = -0.4*(b/L);
a4 = -0.6*(b/L);
a5 = -0.8*(b/L);
a6 = -b/L;

delta_non_p1 = calc_delta_non_parallel(a1);
delta_non_p2 = calc_delta_non_parallel(a2);
delta_non_p3 = calc_delta_non_parallel(a3);
delta_non_p4 = calc_delta_non_parallel(a4);
delta_non_p5 = calc_delta_non_parallel(a5);
delta_non_p6 = calc_delta_non_parallel(a6);

dc_np1 = (12*E*I/(L^3))*(delta_non_p1/F);
dc_np2 = (12*E*I/(L^3))*(delta_non_p2/F);
dc_np3 = (12*E*I/(L^3))*(delta_non_p3/F);
dc_np4 = (12*E*I/(L^3))*(delta_non_p4/F);
dc_np5 = (12*E*I/(L^3))*(delta_non_p5/F);
dc_np6 = (12*E*I/(L^3))*(delta_non_p6/F);

figure
plot(dll, dc_np1);
hold on
plot(dll, dc_np2);
plot(dll, dc_np3);
plot(dll, dc_np4);
plot(dll, dc_np5);
plot(dll, dc_np6);
xlim([0 1]);
ylim([0 1]);
hold off

%% relationship between a tendon’s stiffness and the deflection response to a tip load for both parallel routing and converging routing
r = 0.01;
delta_parallel_tip = (F*L^3)/(3*E*I) - ((r*r*k*L)/(E*I + r*r*k*L))*((F*L^3)/(4*E*I));
delta_converging_tip = (F*L^3)/(3*E*I) - ((b*b*k*L)/(3*E*I + b*b*k*L))*((F*L^3)/(3*E*I));


%% Necessary functions used
function delta_non_parallel = calc_delta_non_parallel(a)
    E = 3500*10^6;
    I = 10^-5;
    L = 0.095;
    z = 0:0.001:0.15;
    F = 1.0;
    b = 0.01;
    delta_non_parallel = (F*z.^3)/(3*E*I) - (((a*z + 3*b).^2).*z)./(3*L*(a*a*L*L + 3*a*b*b + 3*b*b)).*(F*z.^3)/(4*E*I);
end
