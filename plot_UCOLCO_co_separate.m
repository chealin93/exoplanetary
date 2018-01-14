function plot_UCOLCO_ase_co_separate()

close all;

re = load('mr1_re_upto0.9_e8_pi.mat')
unre = load('mr1_unre_upto0.9_e8_pi.mat')

eccen_re = re.eccen;
eccen_unre = abs(unre.eccen);

x_re = 0:0.005:0.7;
x_unre = 0:0.005:0.7;

len1 = eccen_re(1:length(re.LCO_val));
len2 = eccen_unre(1:length(unre.LCO_val));

puco_re = polyfit(eccen_re,re.UCO_val,2)
plco_re = polyfit(len1,re.LCO_val,2)

puco_unre=polyfit(eccen_unre,unre.UCO_val,2)
plco_unre=polyfit(len2,unre.LCO_val,2)

plotpuco_re = polyval(puco_re,x_re);
plotplco_re = polyval(plco_re,x_re);

plotpuco_unre = polyval(puco_unre,x_unre);
plotplco_unre = polyval(plco_unre,x_unre);

figure(1)
    ax=subplot(2,1,1);
    plot(eccen_re,re.UCO_val,'o');hold on;
    plot(x_re,plotpuco_re);
    xlabel('Eccentricity','fontsize',18);
    ylabel('UCO','fontsize',18);
    axis([0 0.7 1 4]);
    bx = subplot(2,1,2);
    plot(len1,re.LCO_val,'o');hold on;
    plot(x_re,plotplco_re)
    xlabel('Eccentricity','fontsize',18)
    ylabel('LCO','fontsize',18)
    axis([0 0.7 1 4]);
figure(2)    
    ax=subplot(2,1,1);
    plot(eccen_unre,unre.UCO_val,'o');hold on;
    plot(x_unre,plotpuco_unre);
    xlabel('Eccentricity','fontsize',18);
    ylabel('UCO','fontsize',18);
    axis([0 0.7 1 4]);
    bx = subplot(2,1,2);
    plot(len2,unre.LCO_val,'o');hold on;
    plot(x_unre,plotplco_unre)
    xlabel('Eccentricity','fontsize',18)
    ylabel('LCO','fontsize',18)
    axis([0 0.7 1 4]);
end


    