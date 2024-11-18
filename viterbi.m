function [path, d] = viterbi(reference, signal, distfun, minback, maxforward)
% [path, d] = viterbi(reference, signal, distfun, windowsize)
%   viterbi alignment of signal to reference
%   reference(:,path) is best match to signal
%   reference is KxD
%   signal is TxD
%   minimizes sum sum_k(distfun(reference(path(k),:), signal(k,:)))
%   subject to constraint that path(k-1) + minback <= path(k) <= path(k-1) + maxforward
%   d is total distance

if (size(reference, 2) > size(reference,1))
    reference = reference';
end

if (size(signal, 2) > size(signal,1))
    signal = signal';
end

%pad reference with Inf on either side to keep it from sticking at a wall
reference = [Inf([1 size(reference,2)]);reference;Inf([1 size(reference,2)])];

K = size(reference, 1);
T = size(signal, 1);

D = zeros(K,T);
P = D;
for i = 1:K
    D(i,1) = distfun(reference(i,:), signal(1,:));
end
tic
for i = 2:T
    for j = 1:K
        i0 = max(j-maxforward,1);
        i1 = max(j-minback,1);
        [mind,ind] = min(D(i0:i1, i-1));
        P(j,i) = ind + i0-1;
        D(j,i) = mind + distfun(reference(j,:), signal(i,:));
    end
%     if (mod(i,ceil(T/20)) == 0)
%         disp([num2str(i) '/' num2str(T)]);
%         toc
%     end
end

path = zeros([T 1]);
[d,path(end)] = min(D(:,T));
for j = T:-1:2
    path(j-1) = P(path(j), j);
end
path = path - 1; %remove extra space from infinity padding

end

