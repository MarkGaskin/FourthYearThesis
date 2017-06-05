%{
Partially Observed Markov Decision Process MATLAB Model

Quantization Section

Developed by: M. Galal, M. Gaskin, I. Harbell, D. Kao

This MATLAB script takes an infinite model and quantizes it
%}

clear all;
close all;
clc;



policy_indicator = 0;       % 0 indicates threshold, 1 indicates all policies



%% Variable Declaration

E= [0 10];              % set Expense cost
R = [0 1];              %set Repair cost
C = 0;

measure_variance = 0;     % error in measurement 
iterative_variance = 1/4;   % breaking on the next step chance
fix_shift = 0.9;     % fixing success
n=100;   %number of belief bins

N=2;    %number of true states 

beta = 0.9;     %discount variable

% d defines the boundaries of the bin, and v assigns a value to each bin
d = [-0.01 1/n:1/n:1];  

%assign values to each bin to be used for C tilde generation
for i = 1:n
    v(i) = (d(i+1)+d(i))/2;
end

%Generate Probabilistic Kernels
[TK,OB] = Generate_Kernels(N, measure_variance, iterative_variance, fix_shift);

%Generate all possible policies, and place in a single matrix
if policy_indicator == 1
    policy = policy_matrix(n);
end

if policy_indicator == 0
    policy = threshold_policy(n);
end


% Calculate number of policies based on policy_indicator
policy_count = n;

if policy_indicator == 1
    policy_count = 2^n;
end

% Declare C tilde matrix
CostMatrix = zeros(N,n);

P_gamma = zeros(n,n);

BK = Belief_Kernel(TK, OB, n, N);

% Generate C tilde
for i = 1:n
    for j = 1:N
        for k=1:N
            CostMatrix(j,i) =  (1-v(i))*E(k) + R(j);
        end
    end
end

J = zeros(n, policy_count);

for k = 1:policy_count
    
    %generate C gamma
    for j = 1:n
        C_gamma(j) = CostMatrix((policy(k,j))+1,j);
    end
    
    %generate P gamma
    for j = 1:n
            P_gamma(j,:) = BK(j,:,policy(k,j)+1);
    end
    
    %Calculate expected cost vector J
    J(:,k) = inv(eye(n)-beta*P_gamma) * transpose(C_gamma);
    
    %Average expected cost vector
    G(1,k) = mean(J(:,k));
    
    % Find the optimal policy
    if k == 1
        Jmin = mean(J(:,k));
        best_policy = policy(k,:);
    end
    
    if mean(J(:,k)) < Jmin
        Jmin = mean(J(:,k));
        best_policy = policy(k,:);
     end
    
end

%Sort the policies
[SortedG, SortedG_index] = sort(G);



