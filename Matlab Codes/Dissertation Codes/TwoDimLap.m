%2D Laplace Equation;
%Ozanski, 14th May

X0=0;
XM=1;
M=8;    %there will be M points in each direction 
            %so (0,.), (x_1,.), ... , (x_M,.) , (XM,.)
        %see jpg file for details
        
h=(XM-X0)/(M+1);  %spacing both horizontally and vertically
%epsilon=1;        

%exact sol:
u = @(x,y) sin(pi*x)*sin(2*pi*y);       

%Defining the coordinates of each triangle (compare with jpg file):
triangles = zeros(2,3,2*(M+1)^2);
% so triangle(:,:,K) =    x1  x2  x3
%                         y1  y2  y3
for j=1:M+1
    for i=1:M+1
        triangles(:,:,2*(j-1)*(M+1)+2*i-1)=[(i-1)*h i*h (i-1)*h    ;                
                                      (j-1)*h   (j-1)*h     j*h];
        triangles(:,:,2*(j-1)*(M+1)+2*i)=[(i-1)*h   i*h       i*h    ;                
                                      j*h   (j-1)*h     j*h];     
    end
end

%Defining the coordinates of each node (see jpg file):
points=zeros(2,(M+2)^2); 
%  so points(:,k) =  [ xk ;
%                      yk];
for i=1:M
    for j=1:M
     points(:,(j-1)*M+i)=[  i*h ;  j*h ];
    end
end

for i=1:M+1
    points(:,M^2+i)=[ (i-1)*h ; 0  ];
    points(:,M^2+M+1+i)=[ XM ; (i-1)*h ];
    points(:,M^2+2*M+2+i)=[ XM-(i-1)*h  ; XM  ];
    points(:,M^2+3*M+3+i)=[ 0 ; XM-(i-1)*h  ];
end

%Deifning the LNODS array (see def in the lecture notes for FEM
%page 39
LNODS=zeros(2*(M+1)^2,3);
for k=1:(2*(M+1)^2)
   for j=1:3
       %so here I consider j-th node of k-th triangle
       for i=1:((M+2)^2)
          if  max(abs(triangles(:,j,k)-points(:,i)))<=h/4  %i.e. j-th node of k-th triangle is i-th point
              LNODS(k,j)=i;
          end
       end
   end
end

%global stiffness matrix:
A=zeros((M+2)^2,(M+2)^2);

%(look at the lecture notes page 39)
for k=1:(2*(M+1)^2)
    CK=K_triangle(triangles(1,1,k),triangles(1,2,k),triangles(1,3,k),...
        triangles(2,1,k),triangles(2,2,k),triangles(2,3,k),epsilon);
    for i=1:3
       for j=1:3
           A(LNODS(k,i),LNODS(k,j))=A(LNODS(k,i),LNODS(k,j))+CK(i,j);
       end        
    end
end

%taking first M^2 rows
%as we test only with phi_i with zero boundary conditions
A=A';
A=A(1:(M^2),1:((M+2)^2));

%Separate the system: values at nodes M^2+1:(M+2)^2 are known from boundary
%cond.
A1=A(:,1:(M^2));
A2=A(:,(M^2+1):((M+2)^2));

U=zeros((M+2)^2,1);
for i=M^2+1:(M+2)^2
U(i)=u(points(1,i),points(2,i));
end
V=U(M^2+1:(M+2)^2);
%finished separating the system.

%Now solve!
U=(A1)\(-A2*V);

%Now plot some results, take without boundary nodes for now:
OUT=zeros(M,M);

for i=1:M
    for j=1:M
OUT(i,j)=U((j-1)*M+i);
    end
end
subplot(1,2,1)
surf(OUT)
title('Numerical Solution')

OUT_exact=zeros(M,M);
for i=1:M
    for j=1:M
OUT_exact(i,j)=u(points(1,(j-1)*M+i),points(2,(j-1)*M+i));
    end
end
subplot(1,2,2)
surf(OUT_exact)
title('Exact Solution')