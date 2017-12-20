%-----------------------------------------------------------------------
% semi_major aixs and period
%-----------------------------------------------------------------------

function [semi_a,period] = kepler3rd(m_t,H)

mul = 1;
adi = 0;
for i =1:length(m_t)
    mul = mul*m_t{i};
    adi = adi+m_t{i};
end

reduced_m = mul/adi;
semi_a = - adi*reduced_m/2/H(1);
period = sqrt(4*pi^2*semi_a^3/adi);
end
