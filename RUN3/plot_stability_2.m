function plot_stability_byangle()
close all;
A = load('mr1_re_upto0.9_e8_pi.mat');

figure(1);
l = 1;
for i = 1:1:length(A.tot_in_dist_cell)
%     xlabel('distance to the center of mass','fontsize',13);
%     ylabel('angle to the center of mass','fontsize',13);
    subplot(7,1,l);
    for j = 1:55
        plot(A.tot_in_dist_cell{i}{j},A.tot_in_ang_cell{i}{j},'ob');hold on;
        plot(A.tot_out_dist_cell{i}{j},A.tot_out_ang_cell{i}{j},'or');
    end
    axis([0.8 3.5 0 180]);
    l = l+1;
end
end

