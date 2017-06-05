function [ policy ] = threshold_policy( n )
% Creates a matrix of all admissible threshold policies given n

k=2^(n-1);
for i=1:n
        x = dec2bin(k,n);
    for j = 1:n
        policy(i,j) = str2num(x(j));
    end
    k=k + 2^(n-1-i);
end
end

