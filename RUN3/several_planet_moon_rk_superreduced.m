
function dy = f(t,y)

global eps;
global m_t;
global tot_num;
global omega;
global R;


n = tot_num;
p = 3;
phi = omega*t;

X=y(1:3);
V=y(4:6);
vel_t = [];
acc_t = [];
dy = [];
%--------position of n-body----------

pos_s = [cosd(90+phi) sind(90+phi) 0]*R;pos_p = -pos_s;
pos_c = X(1:3)';
pos_t = {pos_s pos_p pos_c};

%----------velocity of n-body---------
        
vel_c = V(1:3)';
vel_t = [vel_c];

%----acceleration of thirdbody-----
a3 = - m_t{1}*(pos_t{3}-pos_t{1})/(sqrt(sum((pos_t{3}-pos_t{1}).^2))+eps).^p...
    - m_t{2}*(pos_t{3}-pos_t{2})/(sqrt(sum((pos_t{3}-pos_t{2}).^2))+eps).^p;

acc_t = [a3];

%----------making dy---------------

dy = [vel_t acc_t];
dy = dy';
end