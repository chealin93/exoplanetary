%normalized testparticle mass

function plot_testparticle_energy(H_3,m_t)

nor_H_3=H_3/m_t{3};
length_nor_H_3 = 1:length(nor_H_3);
figure('Name','normalized testparticle energy','NumberTitle','off');
plot(length_nor_H_3,nor_H_3);
end

