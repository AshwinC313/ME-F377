E = 3500*10^6;
I = 10^-5;
L = 0.15;
z = 0:0.001:0.15;
F = 1.0;
delta = (F*z.^3)/(3*E*I);
delta1 = delta - (z/L).*((F*z.^3)/(4*E*I));
dc = (12*E*I/(L^3))*(delta/F);
dc1 = (12*E*I/(L^3))*(delta1/F);
dll = z/L;
figure
plot(dll, dc)
hold on
plot(dll, dc1)