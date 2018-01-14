function plot_UCOLCO_ase_co_our()

close all;


re = load('mr1_re_upto0.7_e8_34pi.mat')

eccen_re = re.eccen;

x_re = 0:0.005:re.eccen(end)+0.2;
x_d = 0:0.03:re.eccen(end)+0.2;

length(x_re)
length(x_d)
syms UCO(e) LCO(e)
UCO(e) = -1.04*e^2 + 2.76*e + 2.37;
LCO(e) = -2.07*e^2 + 2.79*e + 2.09;

format long
UCO_list = UCO(x_d)
LCO_list = LCO(x_d)

len1 = eccen_re(1:length(re.LCO_val));

puco_re = polyfit(eccen_re,re.UCO_val,2)
plco_re = polyfit(eccen_re,re.LCO_val,2)


plotpuco_re = polyval(puco_re,x_re);
plotplco_re = polyval(plco_re,x_re);


figure(1)
    ax=subplot(2,1,1);
    plot(eccen_re,re.UCO_val,'o');hold on;
    plot(x_re,plotpuco_re);
    plot(x_d,UCO_list);
    xlabel('Eccentricity','fontsize',18);
    ylabel('UCO','fontsize',18);
    axis([0 0.7 1.5 4.5]);
    
    ax=subplot(2,1,2);
    plot(eccen_re,re.LCO_val,'o');hold on;
    plot(x_re,plotplco_re);
    plot(x_d,LCO_list);
    xlabel('Eccentricity','fontsize',18);
    ylabel('UCO','fontsize',18);
    axis([0 0.7 1 4.5]);
    
end


    