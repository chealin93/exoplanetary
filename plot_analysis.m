function plot_analysis(T,pos,vel,m_t,D_t,tot_E,tot_U,m)

figure('Name','Analysis','NumberTitle','off');
 
%sun's x-position as a function of time
subplot(3,3,1);
    plot(T,pos{1}(:,1)-pos{2}(:,1));
    title('x-pos of sun');

%position of sun as function of time
subplot(3,3,2);
    plot(T,vel{1}(:,1));
    title('P');
    
%momentum of sun as function of time
subplot(3,3,3);
    plot(T,m_t{1}*vel{1}(:,1));
    title('Q');
    
%distance between two body
subplot(3,3,4);
    plot(T,D_t{1});hold on;
    title('distance of primaries');

%velocity
subplot(3,3,5);
    plot(T,vel{1},'g');hold on;
    plot(T,vel{2},'k');hold on;
%     plot(T,vel{index_moon},'r');
    axis([0 max(T)/10 -1000 1000]);
    title('velocity of S&E&M');

dT = diff(T);
eta(m) = max(dT)/min(dT)
eta_t  = dT/min(dT);

%eta as function of time
subplot(3,3,6);
    plot(T(1:end-1),eta_t);hold on;
    title('Eta as function of t');

%energy and error

H = tot_E + tot_U;

subplot(3,3,7);
    plot(T,tot_E);hold on
    title('kinetic');
subplot(3,3,8);

    plot(T,tot_U);hold on
    title('potential');
subplot(3,3,9);
    plot(T,1-H/mean(H));hold on;
    title('1-tot/mean(tot)');