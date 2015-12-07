function CK = K_triangle( x1,x2,x3,y1,y2,y3,epsilon )
%Computing element stiffness matrix
%Ozanski, 14th May

J=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
J=abs(J);
A123=J/2;

x12=x1-x2;  x21=-x12;
x23=x2-x3;  x32=-x23;
x13=x1-x3;  x31=-x13;
y12=y1-y2;  y21=-y12;
y23=y2-y3;  y32=-y23;
y13=y1-y3;  y31=-y13;
r12=[ x12; y12];  
r23=[ x23; y23];
r13=[ x13; y13];

%the part coming from Laplacian:
A=[ r23'*r23   -r23'*r13   r23'*r12;
    -r23'*r13  r13'*r13    -r12'*r13;
    r23'*r12   -r12'*r13   r12'*r12  ]/(4*A123);

%the non-symmetric part:
p1e=-1;    %=(d phi_1)/(d eta) 
p2e=0;    %=(d phi_2)/(d eta) 
p3e=1;    %=(d phi_3)/(d eta) 
p1x=-1;    %=(d phi_1)/(d xi)
p2x=1;     %=(d phi_2)/(d xi)
p3x=0;     %=(d phi_3)/(d xi)

%see the other jpg file for this:
B=(1/6)*[ y31*p1x-y21*p1e-2*x31*p1x+2*x21*p1e  y31*p1x-y21*p1e-2*x31*p1x+2*x21*p1e  y31*p1x-y21*p1e-2*x31*p1x+2*x21*p1e ;
          y31*p2x-y21*p2e-2*x31*p2x+2*x21*p2e  y31*p2x-y21*p2e-2*x31*p2x+2*x21*p2e  y31*p2x-y21*p2e-2*x31*p2x+2*x21*p2e ;
          y31*p3x-y21*p3e-2*x31*p3x+2*x21*p3e  y31*p3x-y21*p3e-2*x31*p3x+2*x21*p3e  y31*p3x-y21*p3e-2*x31*p3x+2*x21*p3e ];

CK=epsilon*A+B;
end

