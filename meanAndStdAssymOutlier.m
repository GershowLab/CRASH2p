function [u,s] = meanAndStdAssymOutlier(x)

x = x(isfinite(x));

valid = true(size(x));
s = std(x);
u = mean(x);
for j = 1:10
    y = x;
    y(~valid) = u+s;
    u = mean(y);
    s = std(y);
    valid = x < u + s;
end