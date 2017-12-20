%----------------------------------------------------
% plotting eccentricity and energy by each orbit
%----------------------------------------------------

function plot_e_E(orbit_list,e_a_list,ang_a_list,ang_n_list,ang_e_list,e_n_list,del_e_list,total_a_list,total_n_list,energy_error_list)

    figure('Name','discrepancy of analytic and numerical','NumberTitle','off');
    subplot(2,2,1);
        plot(orbit_list,e_a_list,'o');hold on;
        plot(orbit_list,e_n_list,'.','markersize',8);
        
    subplot(2,2,2);
        plot(orbit_list,del_e_list);
    
%     subplot(3,2,3)
%         plot(orbit_list,ang_a_list,'o');hold on;
%         plot(orbit_list,ang_n_list,'.');
%     
%     subplot(3,2,4);
%         plot(orbit_list,ang_e_list);
%         
    subplot(2,2,3);
        plot(orbit_list,total_a_list,'o');hold on;
        plot(orbit_list,total_n_list,'.')
        
    subplot(2,2,4);
        plot(orbit_list,energy_error_list);
end



