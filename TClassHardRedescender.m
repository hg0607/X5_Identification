function y = TClassHardRedescender(x, th)
    if nargin < 2
        th = 3;
    end
    n = length(x);
    y = x;
    for i=1:n
       if abs(x(i)) > th 
           y(i) = 0;
       else
           y(i) = 1;
       end
    end
end
