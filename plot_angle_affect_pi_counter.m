function plot_angle_affect()
close all;
A = load('counter_mr1_re_upto0.9_e4_pi.mat');
B = load('counter_mr1_unre_upto0.9_e4_pi.mat');

distance_range = 0.8:0.05:3.5;
    
inlist0 = []; inlist45 = []; inlist90 = []; inlist135=[]; inlist180=[];

for i = 1:length(A.LCO_val)
    a=0;b=0;c=0;d=0;e=0;
    al = A.LCO_val(i);
    au = A.UCO_val(i);
    
    
    LCO_idx = find(distance_range == al);
    UCO_idx = find(distance_range == au);
    
    
    for j = LCO_idx:1:UCO_idx
        tmp_list = [0 45 90 135 180];
        q=setdiff(tmp_list,A.tot_in_ang_cell{i}{j});
        if length(q) <= 2
            for m = 1:length(q)
                if q(m) == 0
                    a = a +1;
                end
                if q(m) == 45
                    b = b +1;
                end
                if q(m) == 90
                    c = c +1;
                end
                if q(m) == 135
                    d = d +1;
                end
                if q(m) == 180
                    e = e +1;
                end
            end
        end
    end
    inlist0 = [inlist0 a];
    inlist45 = [inlist45 b];
    inlist90 = [inlist90 c];
    inlist135 = [inlist135 d];
    inlist180 = [inlist180 e];        
    
end

inlist0_1 = []; inlist45_1 = []; inlist90_1 = []; inlist135_1=[]; inlist180_1=[];

for i = 1:length(B.LCO_val)
    a=0;b=0;c=0;d=0;e=0;
    al = B.LCO_val(i);
    au = B.UCO_val(i);
    LCO_idx = find(distance_range == al);
    UCO_idx = find(distance_range == au);
    
    for j = LCO_idx:1:UCO_idx
        tmp_list = [0 45 90 135 180];
        q=setdiff(tmp_list,B.tot_in_ang_cell{i}{j});
        if length(q) <= 2
            for m = 1:length(q)
                if q(m) == 0
                    a = a +1;
                end
                if q(m) == 45
                    b = b +1;
                end
                if q(m) == 90
                    c = c +1;
                end
                if q(m) == 135
                    d = d +1;
                end
                if q(m) == 180
                    e = e +1;
                end
            end
        end
    end
    inlist0_1 = [inlist0_1 a];
    inlist45_1 = [inlist45_1 b];
    inlist90_1 = [inlist90_1 c];
    inlist135_1 = [inlist135_1 d];
    inlist180_1 = [inlist180_1 e];        
    
    
end
    
outlist0 = []; outlist45 = []; outlist90 = []; outlist135 =[]; outlist180 =[];

for i = 1:length(A.LCO_val)
    a=0;b=0;c=0;d=0;e=0;
    al = A.LCO_val(i);
    au = A.UCO_val(i);
    LCO_idx = find(distance_range == al);
    UCO_idx = find(distance_range == au);
    
    for j = LCO_idx:1:UCO_idx
        tmp_list = [0 45 90 135 180];
        q=setdiff(tmp_list,A.tot_out_ang_cell{i}{j});
        if length(q) <= 2
            for m = 1:length(q)
                if q(m) == 0
                    a = a +1;
                end
                if q(m) == 45
                    b = b +1;
                end
                if q(m) == 90
                    c = c +1;
                end
                if q(m) == 135
                    d = d +1;
                end
                if q(m) == 180
                    e = e +1;
                end
            end
        end
    end
    
    outlist0 = [outlist0 a];
    outlist45 = [outlist45 b];
    outlist90 = [outlist90 c];
    outlist135 = [outlist135 d];
    outlist180 = [outlist180 e];
end

outlist0_1 = []; outlist45_1 = []; outlist90_1 = []; outlist135_1=[]; outlist180_1=[];
for i = 1:length(B.LCO_val)
    a=0;b=0;c=0;d=0;e=0;
    al = B.LCO_val(i);
    au = B.UCO_val(i);
    LCO_idx = find(distance_range == al);
    UCO_idx = find(distance_range == au);
    
    for j = LCO_idx:1:UCO_idx
        tmp_list = [0 45 90 135 180];
        q=setdiff(tmp_list,B.tot_out_ang_cell{i}{j});
        if length(q) <= 2
            for m = 1:length(q)
                if q(m) == 0
                    a = a +1;
                end
                if q(m) == 45
                    b = b +1;
                end
                if q(m) == 90
                    c = c +1;
                end
                if q(m) == 135
                    d = d +1;
                end
                if q(m) == 180
                    e = e +1;
                end
            end
        end
    end
    
    outlist0_1 = [outlist0_1 a];
    outlist45_1 = [outlist45_1 b];
    outlist90_1 = [outlist90_1 c];
    outlist135_1 = [outlist135_1 d];
    outlist180_1 = [outlist180_1 e];
end

figure(1);
outlist0_1(2)
for i = 1:7
    x = [0 45 90 135 180];
    h(i) = subplot(7,1,i);
        tmp_in = [inlist0(i) inlist45(i) inlist90(i) inlist135(i) inlist180(i)];
        tmp_in1 = [inlist0_1(i) inlist45_1(i) inlist90_1(i) inlist135_1(i) inlist180_1(i)];
        bar(x,[tmp_in' tmp_in1']);
        ylim([0 10])
        
end 
xlabel('degree of angle','fontsize',13);

figure(2);
for i = 1:7
    x = [0 45 90 135 180];
    k(i) = subplot(7,1,i);
        tmp_in = [outlist0(i) outlist45(i) outlist90(i) outlist135(i) outlist180(i)];
        tmp_in1 = [outlist0_1(i) outlist45_1(i) outlist90_1(i) outlist135_1(i) outlist180_1(i)];
        bar(x,[tmp_in' tmp_in1']);
        ylim([0 10])
        
end 
xlabel('degree of angle','fontsize',13);

% save('mr1_re_bf.mat','inlist0')
% save('mr1_re_bf.mat','inlist45','-append')
% save('mr1_re_bf.mat','inlist90','-append')
% save('mr1_re_bf.mat','inlist135','-append')
% save('mr1_re_bf.mat','outlist0','-append')
% save('mr1_re_bf.mat','outlist45','-append')
% save('mr1_re_bf.mat','outlist90','-append')
% save('mr1_re_bf.mat','outlist135','-append')
end
