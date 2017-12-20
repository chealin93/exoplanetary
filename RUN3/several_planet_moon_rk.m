
function dy = f(t,y)

global eps;
global m_t;
global tot_num

n = tot_num;
p = 3;
X=y(1:3*n);

pos_t = {};
j=1;
for i = 1:n
    temp_pos=X(j:j+2);
    pos_t{i} = temp_pos;
    j = j+3;
end

acc_t = [];
k = 1;
for i = 1:length(pos_t)
    exc_pos_t = pos_t;
    exc_m_t = m_t;
    exc_pos_t(i) = [];            %carefully deal with() and {}, they will give different results
    exc_m_t(i) = [];              %https://kr.mathworks.com/help/matlab/matlab_prog/delete-data-from-a-cell-array.html
    acc = 0;
    
    for j = 1: length(exc_pos_t)
        temp_acc = - exc_m_t{j}*(pos_t{i}-exc_pos_t{j})/(sqrt(sum((pos_t{i}-exc_pos_t{j}).^2))+eps).^p;
        acc = acc + temp_acc;
    end
    acc_t(k:k+2) = acc;
    k = k+3;
end

%making dy 
dy = [];
for i = 1:3*n
    dy = [dy y(3*n+i)];
end

for i = 1:3*n
    dy = [dy acc_t(i)];
end

dy = dy';
end