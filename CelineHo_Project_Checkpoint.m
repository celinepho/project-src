% Example of how to checkpoint and restart program using MATLAB. To run this example, just type
% checkpointrestart multiple times.
% To start from begining delete checkpoint file(checkfile.mat)

checkfile='checkfile.mat';
restart = false;
 
%% Set flag to restart program
if(exist(checkfile) == 2)
    restart = true;
else
%% Initialization code
N=11; %set 
dx=1/N;
A_diag=eye(N-1)*(-4/dx^2) %diagonol entries of diag block
A_diag=A_diag+diag(ones(N-2,1),  1)/dx^2; %40000 and 10000 on off diag
A_diag=A_diag+diag(ones(N-2,1), -1)/dx^2; %for diag blocks
A_off=eye(N-1)/dx^2; %construct off diag blocks
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
%Next step is to put blocks in place in big matrix
A=sparse((N-1)*(N-1),(N-1)*(N-1));%Make sparse matrix because matlab can do less calculations-no need to calculate for zeros
for i=1:N-1
    A((i-1)*(N-1)+1:(i-1)*(N-1)+(N-1), (i-1)*(N-1)+1:(i-1)*(N-1)+(N-1))=A_diag; %range of the rows corresponding to ith block do same on right side for columns
end
for i=1:N-1 %for off diag entries, loop again
end
for i=2:N-1 %Starts at 2 since we have 1 less off diag block
     A((i-2)*(N-1)+1:(i-2)*(N-1)+(N-1),(i-1)*(N-1)+1:(i-1)*(N-1)+(N-1))=A_off;
     A((i-1)*(N-1)+1:(i-1)*(N-1)+(N-1),(i-2)*(N-1)+1:(i-2)*(N-1)+(N-1))=A_off; %fill j entry (reverse line 19 i and j)
end

F = zeros(100); %Set f to be zero everywhere to test it
N=100; %number of i and j terms

% F is the 100x1 array
% F is filled up with F(1,1),F(2,1),...F(10,1),F(1,2),F(2,2,),...F(10,2),F(1,3),....,F(9,10),F(10,10)
% F(k) = F(i,j) to populate F where k = i + (j - 1) * 10

for j = 1:N
     y = -pi + 2 / (N-1) * (j - 1) * pi;         
    for i = 1:N
         k = i + (j - 1) * 100;
         x = -pi + 2 / (N-1) * (i - 1) * pi;
         F(k) = sin( pi * (x + pi) / (2 * pi)) * cos( pi / 2 * (2 * (y + pi) / (2 * pi) + 1));
     end
end

%Applying Boundary Conditions
for j = 1 
    x = -pi + 2 / (N-1) * (i - 1) * pi;
    u(i,j)=-pi*(2*pi)^2+(x+pi)*((-2*pi)+(2*pi^2));
   for i = 1 
       y = -pi + 2 / (N-1) * (j - 1) * pi;
     u(i,j)=y*(pi-y)^2;
     for i=10 
         y = -pi + 2 / (N-1) * (j - 1) * pi;
         u(i,j)=(pi-y)^2*cos(y);
           for j=10
              y = -pi + 2 / (N-1) * (j - 1) * pi;
             syms y(u)
             dy=diff(y(u),u);
             dy=0;
         end
     end
   end
end

lambda=1/2;
u=-(A+lambda) \ F; %minus sign to move f on rhs
%u=reshape(u,N-1,N-1) %reshape to make u a 99x99 double
surf(u)
     
%% Increment the main counter
    i=i+1
end