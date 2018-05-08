%Celine Ho
%5/9/2018
%MECE 5397 - Final Project 

clear all;clc

%Construct matrix A
N = 11; %set size of matrix A
dx = 1/N; %distance between each node in x-space
A_diag = eye (N-1) * (-4 / dx^2) %diagonal entries of diagonal block
A_diag = A_diag + diag (ones (N-2,1) , 1) / dx^2; %values of 400 and 100 on off-diag
A_diag = A_diag + diag (ones (N-2,1) , -1) / dx^2; %give values to all diagonal blocks
A_offdiag = eye (N-1) / dx^2; %construct off-diagonal blocks

%Put blocks in place in big matrix
A = sparse ((N-1) * (N-1) , (N-1) * (N-1));%Sparse matrix so matlab can do less calculations-no need to calculate for zeros
for i = 1:N-1
    A((i-1) * (N-1) + 1 : (i-1) * (N-1) + (N-1) , (i-1) * (N-1) + 1 : (i-1) * (N-1) + (N-1)) = A_diag; 
    %range of rows corresponding to ith block do same on right side for columns
end
for i = 1:N-1 %for off-diagonal entries, loop again
end
for i = 2:N-1 %Starts at 2 since we have 1 less off diag block
    %fill i entries
     A((i-2) * (N-1) + 1 : (i-2) * (N-1) + (N-1) , (i-1) * (N-1) + 1 : (i-1) * (N-1) + (N-1)) = A_offdiag;
     %fill j entries
     A((i-1) * (N-1) + 1 : (i-1) * (N-1) + (N-1) , (i-2) * (N-1) + 1 : (i-2) * (N-1) + (N-1)) = A_offdiag; 
end


%Construct matrix F
F = zeros(100); %Set f to be 100x100 array with zeros everywhere to test it
N=100; %number of i and j terms

% F is filled with F(1,1),F(2,1),...F(10,1),F(1,2),F(2,2,),...F(10,2),F(1,3),....,F(9,10),F(10,10)

for j = 1:N
     y = -pi + 2 / (N-1) * (j-1) * pi; %how ‘j’ relates to ‘y’ location       
    for i = 1:N
         k = i + (j-1) * 100; % F(k)=F(i,j) to populate F where k = i + (j - 1) * 10
         x = -pi + 2 / (N-1) * (i-1) * pi; %how ‘i’ relates to ‘x’ location
         F(k) = sin( pi * (x + pi) / (2 * pi)) * cos( pi / 2 * (2 * (y + pi) / (2 * pi) + 1));
     end
end

%Apply Boundary Conditions from problem
for j = 1 
    x = -pi + 2 / (N-1) * (i-1) * pi; %how ‘i’ relates to ‘x’ location
    u(i,j) = -pi * (2*pi)^2 + (x+pi) * ((-2*pi) + (2*pi^2)); %given Dirichlet B.C.
   for i = 1 
       y = -pi + 2 / (N-1) * (j-1) * pi; %how ‘j’ relates to ‘y’ location
     u(i,j) = y * (pi - y)^2; %given Dirichlet B.C.
     for i = 10 
         y = -pi + 2 / (N-1) * (j-1) * pi; %how ‘j’ relates to ‘y’ location
         u(i,j) = (pi - y)^2 * cos(y); %given Dirichlet B.C.
     end
   end
end

m=100; %introduce non-existent grid point (m) for Neumann B.C.

for j = 10
    y = -pi + 2 / (m-1) * (j-1) * pi;
    syms y(u)
    dy=diff (y(u) , u); %given Neumann B.C.
    dy=0; 
end

lambda=1/2;
u=-(A+lambda) \ F; %minus sign to move f onto righthand side
surf(u)
title('Solution to 2D Helmholtz Equation') 