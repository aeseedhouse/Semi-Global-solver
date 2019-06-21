function out = powers(x,m)
% Return matrix of powers of x up to x^m
% 
% inputs:
% x - Nx1 column vector of numbers to take powers of
% m - maximum power
% 
% outputs:
% out - matrix of powers of x
% 
% out =  [1 x(1) x(1)^2 ... x(1)^m
%         1 x(2) x(2)^2 ... x(2)^m
%         ...
%         1 x(N) x(N)^2 ... x(N)^m]

    if size(x,2) > 1
        error('x must be column vector')
    end
    
    out = zeros(length(x),m+1);
    out(:,1) = ones(length(x),1);
    if m > 0
        out(:,2) = x;
        for n=3:m+1
            out(:,n) = out(:,n-1).*x;
        end
    end
end