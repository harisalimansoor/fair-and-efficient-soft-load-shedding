# fair-and-efficient-soft-load-shedding

This is code for the paper "A parametric function based soft load shedding for fair and efficient allocation of electricity to large scale utility customers" . we used a barrier algorithm (a class of interrior point method) and analytical hessian to find an efficient algorithm to solve convex optimization problem with equality and inequality constraitns. The main file contains all the code for results and plots, the following variables should be initialized properly

n=1000; %number of iterations
T1=100; %barrier parameter t
eta=0; %small number for ill to make 1/0 feasible in f, gradient and hessian
stop=10e-4; % stop if variables not change in the last following iterations
beta=0.5; %line search parameter
zeta=0.01; %line search parameter
demand=0.9; % total demand 0.9 mean 90% of total avaialble
no=30; % number of times each experiemtn is performed to computer results

a=[0 0.5 2 10 100 1000 10000];  %alpha=[0-inf] fairness parameter
N=[10,100,1000,10000,100000,1000000]; % number of variables to test
