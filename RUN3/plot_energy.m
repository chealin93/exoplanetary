%----------------------------------------------------
% plotting all body's kinetic and potential energy
%----------------------------------------------------

function plot_energy(T,all_E,all_U)

figure('Name','primary energy','NumberTItle','off');
    subplot(2,3,1);
    plot(T,all_E{1});hold on;
    title('primary kinetic');
    
    subplot(2,3,2);
    plot(T,all_U{1});hold on;
    title('primary potential');
    
    subplot(2,3,3);
    plot(T,all_E{1}+all_U{1});
    title('primary total energy');
    
    subplot(2,3,4);
    plot(T,all_E{2});hold on;
    title('secondary star kinetic');
    
    subplot(2,3,5);
    plot(T,all_E{2});hold on;
    title('secondary potential');
    
    subplot(2,3,6);
    plot(T,all_E{2}+all_U{1});hold on;
    title('secondaty total energy');
end
