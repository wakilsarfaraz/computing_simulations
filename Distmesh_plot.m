addpath distmesh



xmax = 5;
ze = 1;
zm = 4;

%sample_nodes(cx,cy,ze,zm,xmax,density)

[p t] = sample_nodes(2.5,2.5,ze,zm,xmax,0.5);
%simplot(p,t)
%title('Mesh Plot','fontsize',16)

 %  
%   fd=inline('sqrt(sum(p.^2,2))-1','p');
%  [p,t]=distmesh2d(fd,@huniform,0.2,[-1,-1;1,1],[]);
 
 
%  
%  [xx,yy]=meshgrid(-1.1:0.1:1.1,-1.1:0.1:1.1); % Generate grid
%  dd=sqrt(xx.^2+yy.^2)-1; % d(x,y) at grid points >> 
%  [p,t]=distmesh2d(@dmatrix,@huniform,0.2,[-1,-1;1,1],[],xx,yy,dd);
 
 
 
%   fd=inline('-0.56+abs(0.7-sqrt(sum(p.^2,2)))');
%   [p,t]=distmesh2d(fd,@huniform,0.1,[-1,-1;1,1],[]);
%   
%  fd=inline('ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.3))','p');
%  pfix=[-1,-1; -1,1;1,-1;1,1];
%  [p,t]=distmesh2d(fd,@huniform,0.15,[-1,-1;1,1],pfix);

%  fh=inline('min(4*sqrt(sum(p.^2,2))-1,2)','p');
%   [p,t]=distmesh2d(fd,fh,0.05,[-1,-1;1,1],pfix);


%  fd=inline('-0.56+abs(0.9-sqrt(sum(p.^2,2)))');
%  [p,t]=distmesh2d(fd,@huniform,0.1,[-1,-1;1,1],[]);



x = p(:,1);
y = p(:,2);

NNODES = length(x);
NTRI = size(t,1);
LNODES = t;

figure(1)
hold
for i=1:NTRI
    for j=1:3
        trix(j)=x(LNODES(i,j));
        triy(j)=y(LNODES(i,j));
    end
    trix(4)=trix(1);
    triy(4)=triy(1);
    plot(trix,triy)
    title('Distmesh triangulation','fontsize',16)
    xlabel('x','fontsize',14)
    ylabel('y','fontsize',14)
end
axis equal tight
 hold off
 
