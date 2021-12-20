function obj = calAUC( x,datFeat1,datFeat2,lambda )
% Calculate AUC
%   - x: N*M,N is pop size£¬M is dim.
%   - datFeat1: S1*M, Positive
%   - datFeat2: S2*M, Negative
%     tic;

    if nargin < 4
        la = 2.^([-10:1:10]);
        lambda = la(2);
    end
    x1 = datFeat1*x'; %S1*N
    x2 = datFeat2*x'; %S2*N
    N = size(x,1);
    obj = zeros(N,1);
    for i = 1:N
        tmp1 = x1(:,i);
        tmp2 = x2(:,i);
        sum1 = 0;
        for j = 1:size(tmp1,1)
            sum1 = sum1 + sum(tmp1(j) <= tmp2);
        end
        obj(i) = sum1/(size(tmp1,1)*size(tmp2,1));
        obj(i) = obj(i) + 0.5*lambda*norm(x(i,:),2);
    end
%     toc;
end