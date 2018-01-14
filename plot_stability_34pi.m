function plot_stability_byangle_pi()
close all;
A = load('mr1_re_upto0.9_e8_pi.mat');
A.eccen
for i = 1:4:length(A.LCO_val)
    xlabel('Distance to the barycetner[AU]','fontsize',18);
    ylabel('Theta[degree]','fontsize',18);
    axis([0.8 3.5 0 180]);
    figure(i);
    lv = A.LCO_val(i);
    uv = A.UCO_val(i);
    h = area([0.8 lv],[180 180],'Linestyle','none');hold on;
    h.FaceColor = 'r';
    h.FaceAlpha = 0.2;
    j = area([uv 3.5],[180 180],'Linestyle','none');
    j.FaceColor = 'b';
    j.FaceAlpha = 0.2;
    
    for j = 1:55
        plot(A.tot_in_dist_cell{i}{j},A.tot_in_ang_cell{i}{j},'ob');hold on;
        plot(A.tot_out_dist_cell{i}{j},A.tot_out_ang_cell{i}{j},'or');
        
        plot([lv,lv],[0,180],'r','Linewidth',2);
        plot([uv,uv],[0,180],'b','Linewidth',2);
    end
end

    figure(24);
    xlabel('Distance to the barycetner[AU]','fontsize',18);
    ylabel('Theta[degree]','fontsize',18);
    axis([0.8 3.5 0 135]);
    lv = A.LCO_val(24);
    uv = A.UCO_val(24);
    h = area([0.8 lv],[135,135],'Linestyle','none');hold on;
    h.FaceColor = 'r';
    h.FaceAlpha = 0.2;
    j = area([uv 3.5],[135 135],'Linestyle','none');
    j.FaceColor = 'b';
    j.FaceAlpha = 0.2;
       for j = 1:55
        plot(A.tot_in_dist_cell{24}{j},A.tot_in_ang_cell{24}{j},'ob');hold on;
        plot(A.tot_out_dist_cell{24}{j},A.tot_out_ang_cell{24}{j},'or');
        
        plot([lv,lv],[0,135],'r','Linewidth',2);
        plot([uv,uv],[0,135],'b','Linewidth',2);
       end
end