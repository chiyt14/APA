% -------Demostration of the use of APA.m-------
%    Author: Yingtian Chi, Yiwei Qiu
%    Date: Oct, 10, 2020
clear all;
%% ----------------------------------------------
% ------Part 1.single input, single output------
syms x; % define the symbolic variable representing each input
target_fun = @(x) 1 + x(1) + x(1)^2 + exp(x(1)); % define an inline function.
% Also support functions that are defined in seperated files. Please see Part 3.
[fun_approx, OrderAndCoef, collocationPoint, funEval, FunEvalTimes, BasisCoef, BasisCoef_history, LegendreCoef] = APA(target_fun, 1, 1, 1e-2, 6, 0, [], [], [x]);
% ----------OUTPUT------------
% 1. Approximated result
% Input: fun_approx{1}
% Output: ans = 0.043235*x^4 + 0.17613*x^3 + 1.4997*x^2 + 1.998*x + 2.0
% 2. How many times is the target function evaluated ?
% Input: FunEvalTimes
% Output: FunEvalTimes = 13
% 3. What are the sampling points?
% Input: collocationPoint
% Output: collocationPoint = 0    0.5774   -0.5774    0.7746   -0.7746   -0.8611    0.8611   -0.3400    0.3400   -0.9062   -0.5385    0.9062    0.5385
% Explanation: Each column represents a sampling point.
% 4. What are the function values on the sampling points
% Input: funEval
% Output; funEval = 2.0000    3.6920    1.3174    4.5443    1.2863    1.3031    4.9685    1.4874    2.8605    1.3190    1.3351    5.2022    3.5418
% 5. What is the coefficient of each basis function?
% Input: OrderAndCoef{1}
% Output:ans =
%          0    1.0000    2.0000    3.0000    4.0000
%     2.5085    2.1036    1.0245    0.0705    0.0099
% Explanation: It means that the coefficients of 0-th order basis function, ..., and 4-th order basis function are 2.5085, ..., and 0.0099, respectively.
% 6. How to get the function values of the approximated result on the given points (e.g., on these five points: [-1 -0.5 0 0.5 1] )?
% Input: MyLegendre_HD_new([-1 -0.5 0 0.5 1],BasisCoef{1},BasisCoef_history{1},LegendreCoef)
% Output: ans = 1.3688    1.3566    2.0000    3.3986    5.7170
% Comment: For comparison, the real function values of the target function on these points are: 1.3679, 1.3565, 2, 3.3987, 5.7183

%% --------------------------------------------------
% ------Part 2. multiple inputs, single output------
syms x y; % define the symbolic variable representing each input
target_fun = @(x) 1 + x(1) + x(1)^2 + exp(x(2)); % define an inline function.
% Also support functions that are defined in seperated files. Please see Part 3.
[fun_approx, OrderAndCoef, collocationPoint, funEval, FunEvalTimes, BasisCoef, BasisCoef_history, LegendreCoef] = APA(target_fun, 2, 1, 1e-1, 6, 0, [], [], [x y]);
% ----------OUTPUT------------
% 1. Approximated result
% Input: fun_approx{1}
% Output: ans = 5.8739e-16*x^3 + x^2 - 2.4825e-16*x*y + 1.0*x + 0.17394*y^3 + 0.53663*y^2 + 0.99927*y + 1.9963
% 2. How many times is the target function evaluated ?
% Input: FunEvalTimes
% Output: FunEvalTimes = 29
% 3. What are the sampling points?
% Input: collocationPoint
% Output: collocationPoint =
%   1 to 16 columns
%          0    0.5774   -0.5774         0         0    0.5774   -0.5774    0.5774   -0.5774         0         0    0.7746   -0.7746   -0.8611    0.8611   -0.3400
%          0         0         0    0.5774   -0.5774    0.5774    0.5774   -0.5774   -0.5774    0.7746   -0.7746         0         0         0         0         0
%   17 to 29 columns
%     0.3400    0.7746   -0.7746    0.7746   -0.7746    0.5774   -0.5774    0.5774   -0.5774         0         0         0         0
%          0    0.5774    0.5774   -0.5774   -0.5774    0.7746    0.7746   -0.7746   -0.7746   -0.8611    0.8611   -0.3400    0.3400
% Explanation: Each column represents a sampling point.
% 4. What are the function values on the sampling points
% Input: funEval
% Output; funEval =
%   1 to 16 columns
%     2.0000    2.9107    1.7560    2.7813    1.5614    3.6920    2.5373    2.4721    1.3174    3.1697    1.4609    3.3746    1.8254    1.8804    3.6027    1.7756
%   17 to 29 columns
%     2.4556    4.1559    2.6067    2.9360    1.3868    4.0804    2.9257    2.3716    1.2169    1.4227    3.3658    1.7118    2.4049
% 5. What is the coefficient of each basis function?
% Input: OrderAndCoef{1}
% Output:ans =
%          0    1.0000    2.0000    3.0000         0    1.0000    2.0000         0    1.0000         0
%          0         0         0         0    1.0000    1.0000    1.0000    2.0000    2.0000    3.0000
%     2.5085    1.0000    0.6667    0.0000    1.1036   -0.0000         0    0.3578         0    0.0696
% Explanation: The 1st and 2nd rows stand for the order of the 1st input and 2nd input of the basis function, 
%               and the 3rd row stands for the coefficients of the corresponding basis function.
%           For example, the result means that the coefficients of (0,0)-th order basis function, ...,
%               and (0,3)-th order basis function are 2.5085, ..., and 0.0696, respectively.
% 6. How to get the function values of the approximated result on the given points (e.g., on (0,0), (0,0.5) and (-0.5,0.3), 3 points in total)?
% Input: MyLegendre_HD_new([0 0 -0.5; 0 0.5 0.3],BasisCoef{1},BasisCoef_history{1},LegendreCoef)
% Output: ans = 1.9963    2.6519    2.0991
% Comment: For comparison, the real function values of the target function on these points are: 2 2.6487 2.0999



%% ----------------------------------------------------
% ------Part 3. multiple inputs, multiple output------
% use target_fun_demo.m as target functions
[fun_approx, OrderAndCoef, collocationPoint, funEval, FunEvalTimes, BasisCoef, BasisCoef_history, LegendreCoef] = APA(@target_fun_demo, 2, 2, [1e-1 1e-1], 6, 0, [], [], [x y]);
% ----------OUTPUT------------
% ---1. Approximated result---
% Input: fun_approx{1}
% Output: ans = 1.6025e-16*x^2*y + 9.2519e-17*x^2 + x*y + x + 9.2519e-17*y^2 + 1.0*y + 1.0
% Explanation: Approximated result for the 1st output.
% Input: fun_approx{2}
% Output: ans = 0.17394*x^3 + 0.27616*x^2*y^2 + 0.57966*x^2*y + 0.53663*x^2 + 0.57966*x*y^2 + 1.2167*x*y + 0.99927*x + 0.17394*y^3 + 0.53663*y^2 + 0.99927*y + 0.99265
% Explanation: Approximated result for the 2nd output.
% ---2. How many times is the target function evaluated?---
% Input: FunEvalTimes
% Output: FunEvalTimes = 33
% ---3. What are the sampling points?---
% Input: collocationPoint
% Output: collocationPoint =
%   1 to 19 columns
%          0    0.5774   -0.5774         0         0    0.7746   -0.7746    0.5774   -0.5774    0.5774   -0.5774         0         0    0.7746   -0.7746    0.7746   -0.7746   -0.8611    0.8611
%          0         0         0    0.5774   -0.5774         0         0    0.5774    0.5774   -0.5774   -0.5774    0.7746   -0.7746    0.5774    0.5774   -0.5774   -0.5774         0         0
%   20 to 33 columns
%    -0.3400    0.3400    0.5774   -0.5774    0.5774   -0.5774    0.7746   -0.7746    0.7746   -0.7746         0         0         0         0
%          0         0    0.7746    0.7746   -0.7746   -0.7746    0.7746    0.7746   -0.7746   -0.7746   -0.8611    0.8611   -0.3400    0.3400
% Explanation: Each column represents a sampling point.
% ---4. What are the function values on the sampling points---
% Input: funEval
% Output; funEval =
%   1 to 19 columns
%     1.0000    1.5774    0.4226    1.5774    0.4226    1.7746    0.2254    2.4880    0.6667    0.6667    0.1786    1.7746    0.2254    2.7992    0.3555    0.7500    0.0953    0.1389    1.8611
%     1.0000    1.7813    0.5614    1.7813    0.5614    2.1697    0.4609    3.1731    1.0000    1.0000    0.3152    2.1697    0.4609    3.8649    0.8210    1.2180    0.2587    0.4227    2.3658
%   20 to 33 columns
%     0.6600    1.3400    2.7992    0.7500    0.3555    0.0953    3.1492    0.4000    0.4000    0.0508    0.1389    1.8611    0.6600    1.3400
%     0.7118    1.4049    3.8649    1.2180    0.8210    0.2587    4.7077    1.0000    1.0000    0.2124    0.4227    2.3658    0.7118    1.4049
% Explanation: The 1st and the 2nd row stand for the 1st output and 2nd output respectively.
% 5. ---What is the coefficient of each basis function?---
% Input: OrderAndCoef{1}
% Output:ans =
%          0    1.0000    2.0000         0    1.0000    2.0000         0
%          0         0         0    1.0000    1.0000    1.0000    2.0000
%     1.0000    1.0000    0.0000    1.0000    1.0000    0.0000    0.0000
% Explanation: Corresponding to the 1st output. The meaning of each rows is the same as that in Part 2.
% Input: OrderAndCoef{2}
% Output:ans =
%          0    1.0000    2.0000    3.0000         0    1.0000    2.0000         0    1.0000    2.0000         0
%          0         0         0         0    1.0000    1.0000    1.0000    2.0000    2.0000    2.0000    3.0000
%     1.3811    1.2969    0.4191    0.0696    1.2969    1.2167    0.3864    0.4191    0.3864    0.1227    0.0696
% Explanation: Corresponding to the 2nd output. The meaning of each rows is the same as that in Part 2.
% 6. ---How to get the function values of the approximated result on the given points (e.g., on (0,0), (0,0.5) and (-0.5,0.3), 3 points in total)?---
% Input: MyLegendre_HD_new([0 0 -0.5; 0 0.5 0.3],BasisCoef{1},BasisCoef_history{1},LegendreCoef)
% Output: ans = 1.0000    1.5000    0.6500
% Comments: the approximated values of the 1st output of the target function
% Input: MyLegendre_HD_new([0 0 -0.5; 0 0.5 0.3],BasisCoef{2},BasisCoef_history{2},LegendreCoef)
% Output: ans = 0.9927    1.6482    0.7993
% Comments: the approximated values of the 2nd output of the target function
% Comments: For comparison, the real function values of the target function on these points are: (1,1), (1.5,1.6487), (0.6500, 0.8187)


%% ----------------------------------------------------------------------
% ------Part 4. If part of the sampling points have been evaluated------
syms x y; % define the symbolic variable representing each input
target_fun = @(x) 1 + x(1) + x(1)^2 + exp(x(2)); % define an inline function.
[fun_approx, OrderAndCoef, collocationPoint, funEval, FunEvalTimes, BasisCoef, BasisCoef_history, LegendreCoef] = APA(target_fun, 2, 1, 1e-1, 6, 0, [], [], [x y]);
% For example, we have get the approxiamted result of the target function with error control parameter = 1e-1
% fun_approx{1} = 5.8739e-16*x^3 + x^2 - 2.4825e-16*x*y + 1.0*x + 0.17394*y^3 + 0.53663*y^2 + 0.99927*y + 1.9963
% FunEvalTimes = 29
% Now if we are unsatisfied with the accuracy of the result want to set the error control parameter = 1e-2,
% then the sampling points which have been already evaluated can be reused. The method is as follows.
[fun_approx, OrderAndCoef, collocationPoint, funEval, FunEvalTimes, BasisCoef, BasisCoef_history, LegendreCoef] = APA(target_fun, 2, 1, 1e-2, 6, 29, collocationPoint, funEval, [x y]);
% The result is:
% fun_approx{1} = 5.8739e-16*x^3 + x^2 + 4.5179e-16*x*y^3 + 5.3434e-18*x*y^2 - 5.5833e-16*x*y + 1.0*x + 0.043235*y^4 + 0.17613*y^3 + 0.49966*y^2 + 0.99796*y + 2.0
% The number of new sampling points is:
% FunEvalTimes = 12