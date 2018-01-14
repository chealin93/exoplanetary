%-------------------------------------------------------------------
% energy calculation
%-------------------------------------------------------------------

function [all_E,all_U,tot_E,tot_U] = energyc(vel,m_t,D_t,R_m_t,T)

global tot_num

l=length(T);
n = tot_num;
all_E = {};
all_U = {};

%---------------------kinetic-----------------------
for i = 1:n
    E = [];
    for j = 1:l
        kie = 0;
        temp_kie = 1/2*m_t{i}*sum(vel{i}(j,:).^2);
        kie = kie + temp_kie;
        E(j) = kie;
    end
    
    E = E';
    all_E{i} = E;
end

%---------------------potential----------------------
for i = 1:length(D_t)
    U = [];
    for j = 1:l
        poe = 0;
        temp_poe = -R_m_t{i}/D_t{i}(j);
        poe = poe + temp_poe;
        U(j) = poe;
    end
    
    U = U';
    all_U{i} = U;
end

%---total----
tot_E = zeros(size(E));

for i = 1:n
    tot_E = tot_E + all_E{i};
end

tot_U = zeros(size(U));
for i = 1:length(D_t)
    tot_U = tot_U + all_U{i};
end

all_E = all_E;
all_U = all_U;
tot_E = tot_E;
tot_U = tot_U;
end

