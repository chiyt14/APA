function y = target_fun_demo(x)
%    Author: Yingtian Chi, Yiwei Qiu
%    Date: Oct, 10, 2020
%   a 2-outputs-2-inputs target function for demonstrating the use of APA.m
y(1) = 1 + x(1) + x(2) + x(1) * x(2);
y(2) = exp(x(1))*exp(x(2));
end

