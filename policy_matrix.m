function [policy] = policy_matrix(n)
% Creates a matrix of all admissible policies given the number of 
% belief bins n

for i=1:(2^n)
    x = dec2bin(i-1,n);
 
    for j = 1:n
        policy(i,j) = str2num(x(j));
    end
end

end

