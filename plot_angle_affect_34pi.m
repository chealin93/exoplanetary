function plot_angle_affect()
close all;
A = load('mr1_re_upto0.9_e8_pi.mat');
B = load('mr1_unre_upto0.9_e8_pi.mat');

inlist0 = []; inlist45 = []; inlist90 = []; inlist135=[];

for i = 1:length(A.tot_in_ang_cell)
    a=0;b=0;c=0;d=0;
    
    for j = 1:length(A.tot_in_ang_cell{i})
        
        if length(A.tot_in_ang_cell{i}{j}) <= 3
            tmp_list = [0 45 90 135];
            q=setdiff(tmp_list,A.tot_in_ang_cell{i}{j});
            
            for k = length(q)
                if q(k) == 0
                    a = a+1;
                end
                if q(k) == 45
                    b = b+1;
                end
                if q(k) == 90
                    c = c+1;
                end
                if q(k) == 135
                    d = d+1;
                end
            end
        else
            break
        end
    end
    inlist0 = [inlist0 a];
    inlist45 = [inlist45 b];
    inlist90 = [inlist90 c];
    inlist135 = [inlist135 d];
end

inlist0_1 = []; inlist45_1 = []; inlist90_1 = []; inlist135_1=[];

for i = 1:length(B.tot_in_ang_cell)
    a=0;b=0;c=0;d=0;
    
    for j = 1:length(B.tot_in_ang_cell{i})
        
        if length(B.tot_in_ang_cell{i}{j}) <= 3
            tmp_list = [0 45 90 135];
            q=setdiff(tmp_list,B.tot_in_ang_cell{i}{j});
            
            for k = length(q)
                if q(k) == 0
                    a = a+1;
                end
                if q(k) == 45
                    b = b+1;
                end
                if q(k) == 90
                    c = c+1;
                end
                if q(k) == 135
                    d = d+1;
                end
            end
        else
            break
        end
    end
    inlist0_1 = [inlist0_1 a];
    inlist45_1 = [inlist45_1 b];
    inlist90_1 = [inlist90_1 c];
    inlist135_1 = [inlist135_1 d];
end

outlist0 = []; outlist45 = []; outlist90 = []; outlist135 =[];

for i = 1:length(A.tot_out_ang_cell)
    a=0;b=0;c=0;d=0;
    
    for j = 1:length(A.tot_out_ang_cell{i})
        
        if length(A.tot_out_ang_cell{i}{j}) <= 3
            tmp_list = [0 45 90 135];
            q=setdiff(tmp_list,A.tot_out_ang_cell{i}{j});
            
            for k = length(q)
                if q(k) == 0
                    a = a+1;
                end
                if q(k) == 45
                    b = b+1;
                end
                if q(k) == 90
                    c = c+1;
                end
                if q(k) == 135
                    d = d+1;
                end
            end
        else
            break
        end
    end
    outlist0 = [outlist0 a];
    outlist45 = [outlist45 b];
    outlist90 = [outlist90 c];
    outlist135 = [outlist135 d];
end
outlist0_1 = []; outlist45_1 = []; outlist90_1 = []; outlist135_1=[];
for i = 1:length(B.tot_out_ang_cell)
    a=0;b=0;c=0;d=0;
    
    for j = 1:length(B.tot_out_ang_cell{i})
        
        if length(B.tot_out_ang_cell{i}{j}) <= 3
            tmp_list = [0 45 90 135];
            q=setdiff(tmp_list,B.tot_out_ang_cell{i}{j});
            
            for k = length(q)
                if q(k) == 0
                    a = a+1;
                end
                if q(k) == 45
                    b = b+1;
                end
                if q(k) == 90
                    c = c+1;
                end
                if q(k) == 135
                    d = d+1;
                end
            end
        else
            break
        end
    end
    outlist0_1 = [outlist0_1 a];
    outlist45_1 = [outlist45_1 b];
    outlist90_1 = [outlist90_1 c];
    outlist135_1 = [outlist135_1 d];
end

figure(1);
outlist0_1(2)
for i = 1:4
    x = [0 45 90 135];
    h(i) = subplot(4,1,i);
        tmp_in = [inlist0(i) inlist45(i) inlist90(i) inlist135(i)];
        tmp_in1 = [inlist0_1(i) inlist45_1(i) inlist90_1(i) inlist135_1(i)];
        bar(x,[tmp_in' tmp_in1']);
        ylim([0 50])
        
end 
xlabel('degree of angle','fontsize',13);

figure(2);
for i = 1:4
    x = [0 45 90 135];
    k(i) = subplot(4,1,i);
        tmp_in = [outlist0(i) outlist45(i) outlist90(i) outlist135(i)];
        tmp_in1 = [outlist0_1(i) outlist45_1(i) outlist90_1(i) outlist135_1(i)];
        bar(x,[tmp_in' tmp_in1']);
        ylim([0 50])
        
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
