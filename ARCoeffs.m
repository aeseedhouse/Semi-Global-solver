function [b,a] = ARCoeffs(alpha)
    b = zeros(1,2^15);    
    b(1) = 1;
    a = b;
    for i=2:length(b)
        b(i) = (alpha/2+i-1)*b(i-1)/i;
    end
end