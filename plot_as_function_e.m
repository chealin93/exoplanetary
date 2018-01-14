function plot_as_function_e()
close all;

re = load('mr1_re_upto0.9.mat');
unre = load('mr1_unre_upto0.9.mat');

eccen_re = re.eccen;
eccen_unre = unre.eccen;

x_re = 0:0.005:0.1*(length(re.UCO_val)-1);
x_unre = 0:0.005:0.1*(length(unre.UCO_val)-1);

re.LCO_val
re.UCO_val
re.eccen
puco_re=polyfit(eccen_re,re.UCO_val,2)
plco_re=polyfit(eccen_re,re.LCO_val,2)

puco_unre=polyfit(eccen_unre,unre.UCO_val,2)
plco_unre=polyfit(eccen_unre,unre.LCO_val,2)

plotpuco_re = polyval(puco_re,x_re);
plotplco_re = polyval(plco_re,x_re);

plotpuco_unre = polyval(puco_unre,x_unre);
plotplco_unre = polyval(plco_unre,x_unre);

figure(1)
    ax=subplot(2,1,1);
    plot(eccen_re,re.UCO_val,'o');hold on;
    plot(x_re,plotpuco_re)
    xlabel('eccentricity','fontsize',13)
    ylabel('UCO','fontsize',13)
    title('restricted problem','fontsize',15)
    
    bx = subplot(2,1,2);
    plot(eccen_re,re.LCO_val,'o');hold on;
    plot(x_re,plotplco_re)
    xlabel('eccentricity','fontsize',13)
    ylabel('LCO','fontsize',13)
    
figure(2)
    subplot(2,1,1);
    
    plot(eccen_unre,unre.UCO_val,'o');hold on;
    plot(x_unre,plotpuco_unre)
    xlabel('eccentricity','fontsize',13)
    ylabel('LCO','fontsize',13)
    title('restricted problem','fontsize',15)
    
    subplot(2,1,2)
    plot(eccen_unre,unre.LCO_val,'o');hold on;
    plot(x_unre,plotplco_unre)
    xlabel('eccentricity','fontsize',13)
    ylabel('LCO','fontsize',13)
    
end


    