% This is an example script to set up and solve a parametric optimization problem 
% using GUROBI, for a range of values of a scalar parameter (gamma).

clear all; clc;
 
nx = 1; % no. of decision variables

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% DEFINING AFFINE MODES (nq is the number of affine moodes)

% nq = 2;
% A = [-1;1];
% B = [0;0];

% nq = 3;
% A = [-1;1;0.5];
% B = [0;0;1];

nq = 4;
A = [-1;1;0;0.5];
B = [0;0;1.5;1]; 
 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% DEFINING THE OPTIMIZATION PROBLEM

nz = 4.*nx+4+3.*nq;
na = 3 + 3.*nq;
gamma = 0:0.01:5;
 
lb = -2;
ub = 3;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% PLOTTING (FOR INTUITIVE PURPOSES)

x = lb:0.01:ub;
y = A*x+B;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% BUILDING GUROBI MODEL
 
for i=1:length(gamma)
 
    model.obj = zeros(nz,1);
    model.obj(end-nx-3:end-nx-1) = [0.5;0.5;-1];
 
    model.lb = [repmat(lb,3,1);-Inf(nz-3.*nx,1)];
    model.ub = [repmat(ub,3,1);Inf(nz-3.*nx,1)];
 
    model.A = sparse(zeros(na,nz));
    model.rhs = zeros(na,1);
    model.sense = repmat('=',na,1);
   
    model.rhs(1:3) = [0;0;gamma(i)];
    model.rhs(4:end) = repmat(B,3,1);
 
    model.A(1:2,1:3) = [-0.5,-0.5,1;-1,1,0];
    model.A(2,end-1) = 1;
    model.A(3,end) = 1;
 
    model.A(4:end,1:3*nx) = blkdiag(-A,-A,-A);
    model.A(4:end,3*nx+1:end-4-nx) = eye(3*nq);

   
    model.genconnorm.resvar = nz;
    model.genconnorm.vars = nz-1;
    model.genconnorm.which = 2;
    
    model.genconmax(1).resvar = 3.*nx+3.*nq+1;
    model.genconmax(1).vars = 3.*nx+(1:nq);
 
    model.genconmax(2).resvar = 3.*nx+3.*nq+2;
    model.genconmax(2).vars = 3.*nx+nq+(1:nq);
 
    model.genconmax(3).resvar = 3.*nx+3.*nq+3;
    model.genconmax(3).vars = 3.*nx+nq+nq+(1:nq);
 
 
    out=gurobi(model);
    J(i) = out.objval;
end

cl
% plot f vs x
% plot(x,max(y));cle

% to plot h1 vs. gamma
% plot(gamma(1:500),J)
