function plot_UCOLCO_ase_co()

close all;

re = load('mr1_re_upto0.7_e8_pi_1.mat')

eccen_UCO = re.eccen_UCO;
eccen_LCO = re.eccen_LCO;

x_re = 0:0.005:0.7;

puco_re = polyfit(eccen_UCO,re.UCO_val,2)
plco_re = polyfit(eccen_LCO,re.LCO_val,2)

plotpuco_re = polyval(puco_re,x_re);
plotplco_re = polyval(plco_re,x_re);

figure(1)
    ax=subplot(1,1,1);
    plot(eccen_UCO,re.UCO_val,'o');hold on;
    plot(x_re,plotpuco_re);
    plot(eccen_LCO,re.LCO_val,'o');hold on;
    plot(x_re,plotplco_re);
    xlabel('Eccentricity','fontsize',18);
    ylabel('UCO','fontsize',18);
    axis([0 0.7 0 5]);
end


    