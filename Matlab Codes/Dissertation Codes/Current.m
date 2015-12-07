%%  This program solves the limiting case of the model equation (Simpler Problem).
%   Finite Element method is coded in Cylindrical polar coordinates
%   to solve for partial pressure of oxygen.
%    Wakil Sarfaraz. 22 May 2014 
clear all;
tic
xmax = 10;
N = 64;%(Number of points on the x, y interval on which the equation is solved.
X = linspace(0,xmax,N+1); %This divides the interval into N equispaced sub intervals.
%X = 0: 1/N :1; This can also be used to create the same X.
[x, y] = meshgrid(X,X); % This creates an (N+1) by (N+1) grid of values ...
                        % of X and assignes them to x and y. (So x and y
                        % are now matrices of size (N+1)x(N+1)
x =x(:);  % This turns the matrix x into a column vector.
y =y(:);  % This turns the matrix y into a column vector.
NNODES = (N+1)^2; % This introduces the number of nodes in the result of triangulation.
NTRI = 2*N^2;  % This is the number of triangles in the mesh.
LNODES = zeros(NTRI, 3); % This creates a zero matrix of the size 'number of triangles by 3'
                         % (since there are 3 vertices for each triangle.

for i = 1:N
    for j = 1: N
        LNODES(i+2*(j-1)*N,1) = i+(j-1)*(N+1);%This numbers each of the first nodes of all the lower triangles in each square.
        LNODES(i+2*(j-1)*N,2) = i+j*(N+1); % This numbers each of the second nodes of all the lower triangles in each square.
        LNODES(i+2*(j-1)*N,3) = (i+1)+(j-1)*(N+1); % This numbers each of the third nodes of all the lower triangles in each square.
        
        LNODES(i+N+2*(j-1)*N,1) = i+1+j*(N+1); %This numbers each of the first nodes of all the upper triangles in each square.
        LNODES(i+N+2*(j-1)*N,2) = (i+1)+(j-1)*(N+1);%This numbers each of the second nodes of all the upper triangles in each square.
        LNODES(i+N+2*(j-1)*N,3) = i+j*(N+1);%This numbers each of the third nodes of all the upper triangles in each square.
                                            % Here for those nodes that
                                            % some triangles share (copy
                                            % and paste can work).
    end
end

SP = sparse(NNODES, NNODES);  %This creates a sparse matrix of the size 'Number of nodes by Number of nodes'(entries usually computed by integration).
LV = zeros(NNODES,1);         % This creates the load vector of zero entries for now which is l(v) part of the weak formulation.

for n = 1: NTRI
    r1 = [x(LNODES(n,1)) y(LNODES(n,1))];% Position vector (x,y)' for the first nodes of all triangles.
    r2 = [x(LNODES(n,2)) y(LNODES(n,2))];% Position vector (x,y)' for the second nodes of all triangles.
    r3 = [x(LNODES(n,3)) y(LNODES(n,3))];% Position vector (x,y)' for the third nodes of all triangles.
    J = [r2(1)-r1(1) r2(2)-r1(2); r3(1)-r1(1) r3(2)-r1(2)]; %This is the jacobian matrix of the mapping.

    detJ(n) = det(J);
   Astiff = ((r2(1)+r3(1))/(2*det(J)))*...
       [(r2(2)-r3(2))^2+(r3(1)-r2(1))^2 ...
       (r2(2)-r3(2))*(r3(2)-r1(2))+(r3(1)-r2(1))*(r1(1)-r3(1))...
       (r2(2)-r3(2))*(r1(2)-r2(2))+(r3(1)-r2(1))*(r2(1)-r1(1)); ...
       (r3(2)-r1(2))*(r2(2)-r3(2))+(r1(1)-r3(1))*(r3(1)-r2(1)) ...
       (r3(2)-r1(2))^2+(r1(1)-r3(1))^2  ...
       (r3(2)-r1(2))*(r1(2)-r2(2))+(r1(1)-r3(1))*(r2(1)-r1(1)); ...
       (r2(2)-r3(2))*(r1(2)-r2(2))+(r3(1)-r2(1))*(r2(1)-r1(1))  ...
       (r3(2)-r1(2))*(r1(2)-r2(2))+(r1(1)-r3(1))*(r2(1)-r1(1))  ... 
       (r1(2)-r2(2))^2+(r2(1)-r1(1))^2];
       for i = 1 : 3 
           for j = 1:3  
               SP(LNODES(n,i),LNODES(n,j)) =SP(LNODES(n,i),LNODES(n,j))...
                   + Astiff(i,j); %The values of Astiff in the SP.
           end
       end
     
end
for i = 1: NNODES
    if (x(i)==xmax || y(i)==xmax)  % This enforces boundary conditions.
        LV(i) = 1;
        SP(i,:) = 0;
        SP(i,i) =1 ;
     elseif (y(i) == 0 && x(i) >= 0 && x(i) <= 1)
        SP(i,:) = 0;
        SP(i,i) = 1;
        LV(i) = 0;
    end
end

U = SP\LV; %Solves the linear system.

trisurf(LNODES,x,y,U)
%shading interp
xlabel('r','fontsize',14) 
xlim([0 xmax])
ylim([0 xmax])
ylabel('z','fontsize',14)
zlabel('q(r,z)','fontsize',14)
title('Finite Element Solution for partial pressure','fontsize',18)


% Computation of the Current using finite different approximation,
DU = zeros((N+1)^2,1);

I =0;
for i = 1: (N+1)^2
    
    if y(i)==0 && x(i)<1   
     DU(i) = (U(i+1)-U(i))/(y(i+1)-y(i));
    I = I+2*x(i)*DU(i);
    elseif y(i)==0 && x(i)==1
      DU(i) = (U(i+1)-U(i))/(y(i+1)-y(i));
    I = I+DU(i);   
    end
    
end
I = (xmax/(N+1))*pi/4*I


% Computation of the current using the Finite element triangles.
Dh = zeros(NTRI,1);
Ih = 0;


for i = 1: NTRI
   
    Dh(i) = ([x(LNODES(i,1))-x(LNODES(i,3)) y(LNODES(i,1))-y(LNODES(i,3))]*...
        [x(LNODES(i,1))-x(LNODES(i,3)) y(LNODES(i,1))-y(LNODES(i,3))]'...
        *(U(LNODES(i,2))-U(LNODES(i,1)))*(U(LNODES(i,1))-U(LNODES(i,2)))...
        +[x(LNODES(i,1))-x(LNODES(i,2)) y(LNODES(i,1))-y(LNODES(i,2))]*...
        [x(LNODES(i,1))-x(LNODES(i,2)) y(LNODES(i,1))-y(LNODES(i,2))]'...
        *(U(LNODES(i,3))-U(LNODES(i,1)))*(U(LNODES(i,1))-U(LNODES(i,3)))...
        +[x(LNODES(i,2))-x(LNODES(i,1)) y(LNODES(i,2))-y(LNODES(i,1))]*...
        [x(LNODES(i,1))-x(LNODES(i,3)) y(LNODES(i,1))-y(LNODES(i,3))]'...
        *((U(LNODES(i,3))-U(LNODES(i,1)))*(U(LNODES(i,1))-U(LNODES(i,2)))...
        +(U(LNODES(i,2))-U(LNODES(i,1)))*(U(LNODES(i,1))-U(LNODES(i,3)))))...
        *(1/6)*(x(LNODES(i,1))+x(LNODES(i,2))+x(LNODES(i,3)));
    
    r1 = [x(LNODES(i,1)) y(LNODES(i,1))];% Position vector (x,y)' for the first nodes of all triangles.
    r2 = [x(LNODES(i,2)) y(LNODES(i,2))];% Position vector (x,y)' for the second nodes of all triangles.
    r3 = [x(LNODES(i,3)) y(LNODES(i,3))];% Position vector (x,y)' for the third nodes of all triangles.
    J = [r2(1)-r1(1) r2(2)-r1(2); r3(1)-r1(1) r3(2)-r1(2)]; %This is the jacobian matrix
                                              
    Ih = Ih+Dh(i)/det(J);                    
      
  
end
Ih = (-pi/2)*Ih

Differ = abs(I-Ih);
TI = abs(1-I);
TIh = abs(1-Ih);
toc

