function plot_test_dvorak()

syms e;
e = 0:0.05:0.7;

UCO_d = 2.37 + 2.76*e - 1.04*e.^2;
LCO_d = 2.09 + 2.79*e - 2.07*e.^2;

UCO_re = 2.3022 + 3.012*e - 2.0278*e.^2;
LCO_re = 1.9914 + 1.1001*e - 1.9747*e.^2;

UCO_unre = 2.2709 + 2.9433*e - 1.7776*e.^2;
LCO_unre = 1.9636 + 1.3195*e - 2.2573*e.^2;


figure(1);
subplot(2,1,1);
plot(e,UCO_d);hold on;
plot(e,UCO_re);
plot(e,UCO_unre);
subplot(2,1,2);
plot(e,LCO_d);hold on;
plot(e,LCO_re);
plot(e,LCO_unre);
end


