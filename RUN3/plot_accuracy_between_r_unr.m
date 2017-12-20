function plot_accuracy_between_r_unr()

A=load('restricted.mat');
A1 = load('restricted1.mat');
B=load('unrestricted.mat');
B1 = load('unrestricted1.mat');
ah = length(A.H3)
bh = length(B.H3)

ah1 = length(A1.H3)
bh1 = length(B1.H3)

a=A.H3(1:bh)-B.H3;
b=A1.H3(1:bh1)-B1.H3;

a1 = 1:length(a);
b1 = 1:length(b);
figure('Name','accuracy between restricted and unrestricted','NumberTitle','off');

ax = subplot(1,1,1)
    plot(a1,a);hold on;
    plot(b1,b);
    xlabel('Number of Integration');
    ylabel('Total energy of Third body');
    xlim([0 500]);

    
end
