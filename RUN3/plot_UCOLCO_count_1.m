function plot_UCOLCO_count()

close all;

re = load('counter_mr1_re_upto0.7_e8_pi_1.mat');
unre = load('counter_mr1_unre_upto0.7_e8_pi_1.mat');

eccen_re_UCO = re.eccen_UCO;
eccen_re_LCO = re.eccen_LCO
eccen_unre_UCO = abs(unre.eccen_UCO);
eccen_unre_LCO = abs(unre.eccen_LCO);

x_re = 0:0.005:0.71;
x_unre = 0:0.005:0.71;

puco_re = polyfit(eccen_re_UCO,re.UCO_val,2)
plco_re = polyfit(eccen_re_LCO,re.LCO_val,2)

puco_unre=polyfit(eccen_unre_UCO,unre.UCO_val,2)
plco_unre=polyfit(eccen_unre_LCO,unre.LCO_val,2)

plotpuco_re = polyval(puco_re,x_re);
plotplco_re = polyval(plco_re,x_re);

plotpuco_unre = polyval(puco_unre,x_unre);
plotplco_unre = polyval(plco_unre,x_unre);

figure(1)
    ax=subplot(2,1,1);
    plot(eccen_re_UCO,re.UCO_val,'o');hold on;
    plot(x_re,plotpuco_re);
    plot(eccen_unre_UCO,unre.UCO_val,'o');
    plot(x_unre,plotpuco_unre);
    axis([0 0.71 0 4])
    xlabel('eccentricity','fontsize',13)
    ylabel('UCO','fontsize',13)
    title('restricted problem','fontsize',15)
    
    bx = subplot(2,1,2);
    plot(eccen_re_LCO,re.LCO_val,'o');hold on;
    plot(x_re,plotplco_re)
    plot(eccen_unre_LCO,unre.LCO_val,'o');hold on;
    plot(x_unre,plotplco_unre);
    xlabel('eccentricity','fontsize',13)
    ylabel('LCO','fontsize',13)
    axis([0 0.71 0 4])
end

