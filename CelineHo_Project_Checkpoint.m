%% checkpoint and restart program

checkfile='checkfile.mat';
restart = false;
 
%% Set flag to restart program
if(exist(checkfile) == 2)
    restart = true;
else
%% Initialization code
N = 11; %set number of nodes which effects size of matrix A
dx = 1/N; %distance between each node in x-space
onesmatrix = ones (N-2,1); %9x1 matrix of ones

A_diag = eye (N-1) * (-4 / dx^2) %initial value (-4) for diag entries of diag block
A_diag = A_diag + diag (onesmatrix , 1) / dx^2; %values for 1st diagonal of diag block
A_diag = A_diag + diag (onesmatrix , -1) / dx^2; %values for diag below main diagonal of diag blocks
end
 
%% Main loop written as a while loop
i = 1;
while (i < 1000)
    %% Save the state and exit
    if( mod(i,10) == 0)
        disp('Checkpointing Program');
        save(checkfile);
        return;
    end
 
    %% Restart the program if the restart flag is set
    if(restart)
        disp('Restarting Program');
        load(checkfile);
        restart = false;
    end
     
%% Main body of the loop
A_offdiag = eye (N-1) / dx^2; %construct identity matrix for off-diag blocks
%Put blocks in place in big matrix
A = sparse ((N-1) * (N-1) , (N-1) * (N-1));%Sparse 100x100 matrix so less calculations needed
%no need to calculate for zeros

for i = 1:N-1
    A((i-1) * (N-1) + 1 : (i-1) * (N-1) + (N-1) , (i-1) * (N-1) + 1 : (i-1) * (N-1) + (N-1)) = A_diag; 
   %range of rows&columns corresponding to ith block
   %fill each row&column value in A_diag block
end
for i = 2:N-1 %Starts at 2 since there is 1 less off diag block
     %fill upper diag blocks
     A((i-2) * (N-1) + 1 : (i-2) * (N-1) + (N-1) , (i-1) * (N-1) + 1 : (i-1) * (N-1) + (N-1)) = A_offdiag;
     %%fill lower diag blocks (switch i&j entries)
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
     
%% Increment the main counter
    i=i+1
end