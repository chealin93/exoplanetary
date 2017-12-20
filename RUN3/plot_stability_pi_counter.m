function plot_stability_pi_counter()
close all;
A = load('counter_mr1_re_upto0.9_e8_pi.mat');

for i = 1:floor(length(A. LCO_val)/5)
    figure(i);
    xlabel('Distance to the barycetner[AU]','fontsize',18);
    ylabel('Theta[degree]','fontsize',18);
    h  = 0;
    for j = i*5-4:i*5
        subplot(5,1,h+1);
        axis([0.8 3.5 0 180]);hold on;
        lv = A.LCO_val(j);
        uv = A.UCO_val(j);
        
        p = area([0.8 lv],[180 180],'Linestyle','none');hold on;
        p.FaceColor = 'r';
        p.FaceAlpha = 0.2;
        q = area([uv 3.5],[180 180],'Linestyle','none');
        q.FaceColor = 'b';
        q.FaceAlpha = 0.2;
    
        for k = 1:55
            plot(A.tot_in_dist_cell{j}{k},A.tot_in_ang_cell{j}{k},'ob');hold on;
            plot(A.tot_out_dist_cell{j}{k},A.tot_out_ang_cell{j}{k},'or');
            
            plot([lv,lv],[0,180],'r','Linewidth',2);
            plot([uv,uv],[0,180],'b','Linewidth',2);
        end
        h = h+1;
    end
end

end