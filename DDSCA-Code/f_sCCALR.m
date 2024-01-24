function [w1, w2, obj] = f_sCCALR(X, Y, Z, paras)


% set parameters
t1 = paras.lambda1;
t2 = paras.lambda2;
s1 = 1;
s2 = 1;

% initialize canonical loadings
n_XVar = size(X,2);
n_YVar = size(Y,2);

w1 = ones(n_XVar, 1)./n_XVar;
w2 = ones(n_YVar, 1)./n_YVar;


% stop criteria
stop_err = 10e-5;


max_iter = 50;
for iter = 1:max_iter
    
    % fix w2, get w1
    res = Y*w2;   
    XX = X'*X;
    XY = X'*res;    
    Wi = sqrt(sum(w1.*w1,2)+eps);
    D1 = diag(1./Wi);
    w1 = ((1+s1)*XX+t1*D1)\XY;
    scale1 = sqrt((X*w1)' * X * w1);
    w1 = w1 ./ scale1;
    Xw1 = X*w1;
    
    % fix w1, get w2
    res = X*w1;
    YY = Y' * Y;
    XY = Y' * res;
    Wi = sqrt(sum(w2.*w2,2)+eps);
    D1 = diag(1./Wi);
    w2 = ((2+s2)*YY+t2*D1)\(XY+Y'*Z);
    scale2 = sqrt((Y*w2)' * Y * w2);
    w2 = w2 ./ scale2;
    Yw2 = Y*w2;
    
    obj(iter) = abs(norm(Yw2-Xw1,2));
    
    if iter > 2 && abs(obj(iter) - obj(iter-1)) < stop_err
        break;
    end
    
end
