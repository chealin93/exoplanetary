%----------------------------------------------------
% plotting eccentricity and energy by each orbit
%----------------------------------------------------

function plot_e_E1(orbit_list,e_a_list,ang_a_list,ang_n_list,ang_e_list,e_n_list,del_e_list,total_a_list,total_n_list,energy_error_list)
    
    figure('Name','discrepancy of analytic and numerical of eccentricity','NumberTitle','off');
        subplot(2,1,1);
            plot(orbit_list,e_a_list,'o');hold on;
            plot(orbit_list,e_n_list,'.','markersize',8);
            xlabel('orbits','fontsize',13);
            ylabel('eccentricity','fontsize',13);
        subplot(2,1,2);
            plot(orbit_list,del_e_list);
            xlabel('orbits','fontsize',13);
            ylabel('discrepancy of eccentricity','fontsize',13);
        
    figure('Name','discrepancy of analytic and numerical of total energy','NumberTitle','off');
        
        subplot(2,1,1);
            format compact
            plot(orbit_list,total_a_list,'o');hold on;
            plot(orbit_list,total_n_list,'.')
            xlabel('orbits','fontsize',13);
            ylabel('total energy','fontsize',13);
        
        subplot(2,1,2);
            plot(orbit_list,energy_error_list);
            xlabel('orbits','fontsize',13);
            ylabel('discrepancy of total energy','fontsize',13);
      
end




