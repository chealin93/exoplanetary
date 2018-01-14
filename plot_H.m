%----------------------------------------------------
% plotting Hamiltonian
%----------------------------------------------------

function plot_H(T,tot_E,tot_U)

figure('Name','entire system"s energy','NumberTitle','off');
    subplot(3,1,1);
    plot(T,tot_E);
    title('total kinetic');
    subplot(3,1,2);
    plot(T,tot_U);
    title('total potential');
    subplot(3,1,3);
    plot(T,tot_E+tot_U);
    title('hamiltonian')
end

