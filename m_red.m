function [m_red] = m_red(m_t)

mul = 1;
adi = 0;
for i = 1:length(m_t)
    mul = mul*m_t{i};
    adi = adi + m_t{i};
end

m_red = mul/adi;
    
    