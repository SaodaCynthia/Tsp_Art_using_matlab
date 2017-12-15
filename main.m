clc;
global magnify 
global neg_color
global pos_color
global conv_limit
global def_resolution
 neg_color=255;
 pos_color= 0;
 conv_limit=.5;
 def_resolution=1;
 magnify=10;

[I]=imread('smile.jpg'); 
size(I)
figure(1)
imshow(I)
title('input image');
[I]=rgb2gray(I);
figure(2)
imshow(I)
title('grayscale version input image');
[I]=imadjust(I,[.3 .8],[]);
figure(3)
imshow(I)
title('increased contrast of grayscale');
input= imresize(I,[300,300]);
title('resized grayscale');
figure(4)
imshow(input)

size(input)
sImage=voronoi_stipple(input);
figure(6)
imshow(sImage)
title('TSP art');
hold on

[resize]=imresize(sImage,[300 300]);
[h,w]=size(resize);
nodes=[];

for i=1:w
    for j=1:h
       pix=impixel(sImage,j,i);
       pixel(j,i)=pix(1);
       if i==1 && j==1
       neg_color=pixel(j,i);
       end
       if pixel(j,i)>=0 && pixel(j,i)< 30
       nodes=[nodes;[j,i]];
       end
    end
end
fprintf('Read in %d nodes', length(nodes));
[lanes,P]=connect_dots_with_lines(nodes);
x = nodes(P,1);
x = [x;x(1)];
y = nodes(P,2);
y = [y;y(1)];
plot(x,y,'k')
axis equal


function [sImage]= voronoi_stipple (im) 
global magnify 
global neg_color
global pos_color
global conv_limit
global def_resolution
 neg_color=255;
 pos_color= 0;
 conv_limit=.5;
 def_resolution=1;
 magnify=8;

[height,width]= size (im);
imsize=[height,width];


n_cells= (floor(hypot(height,(width/3))))* magnify;

fprintf('Creating %d stipples with %d conevergence point',n_cells,conv_limit);

centroids= zeros(2,n_cells);
for i=1:n_cells
  for j=1:2
    a=randi(height);
    
    b=randi(width);
    if j==1
    centroids(j,i)= a;
    else
    centroids(j,i)= b;
    end
  end
end
rho=zeros(width,height);
pixel=zeros(width,height);

for y=1:width
    for x=1:height
        pix = impixel(im,x,y);
        pixel(y,x)=pix(1);
        rho(y,x)= 1-pix(1)/255;
        
    end
end

[pixel]=clr_image(height,width,pixel);
fImage=draw_points(centroids,pixel); 
% imshow(fImage);
% imwrite(fImage,'vstipple.png');

newcentroids_sum=zeros(3,n_cells);
iter=1; resolution= def_resolution;
res_inv=1/def_resolution;

while(iter<300)
    [newcentroids_sum]=zero_sums(newcentroids_sum);
    [newcentroids_sum]=sumofregions(centroids, newcentroids_sum, rho,res_inv ,height,width);
    [centroidal_del,centroids,newcentroids_sum]= compute_centroids(centroids, newcentroids_sum, imsize);
    fprintf('%d \tdifference %f', iter, centroidal_del);
    clr_image(height,width,pixel);
    sImage=draw_points(centroids,pixel);
    size(sImage);
%     figure(4+iter)
%     imshow(sImage);
  
    if centroidal_del==0
        resolution =resolution*2;
    elseif centroidal_del< conv_limit*resolution
        break;
    end
    iter=iter+1;
end

% [tImage]=draw_pointswith_magnify(centroids,height,width);
figure(5)
imshow(sImage);
end

function [centroidal_del,centroids,newcentroids_sum]= compute_centroids (centroids,newcentroids_sum,imsize)

centroidal_del=0;
[~,c]=size(centroids);
for i= 1:c
    if ~newcentroids_sum(3,i)
        centroids(1,i)=randi(imsize(1));
        centroids(2,i)=randi(imsize(2));
        
    else
        newcentroids_sum(1,i)=newcentroids_sum(1,i)/newcentroids_sum(3,i);
        newcentroids_sum(2,i)=newcentroids_sum(2,i)/newcentroids_sum(3,i);
        centroidal_del=centroidal_del+ hypot((newcentroids_sum(1,i)-centroids(1,i)),(newcentroids_sum(2,i)-centroids(2,i)));
        centroids(1,i)=newcentroids_sum(1,i);
        centroids(2,i)=newcentroids_sum(2,i);
    end
end
end

function [newcentroids_sum]= sumofregions(centroids,newcentroids_sum,rho,resolution_step,height,width)
global points; 

h_range=resolution_step/2:resolution_step:height;
w_range=resolution_step/2:resolution_step:width;
points= cartesian(h_range,w_range);
[r,c]=size(points);
centro=centroids';
nearest_neighbours_indices=zeros(1,r);
nearest_neighbours_indices=knnsearch(centro,points,'k',1);

for i=1:length(points)
 %nearest_neighbours_indices(i)=knnsearch(centro,points(i,:),'k',1);
 pt=points(i,:);
 xa=pt(1);ya=pt(2);
 x=ceil(xa);y=ceil(ya);
 r=rho(x,y);
 nn_index=nearest_neighbours_indices(i);
 newcentroids_sum(1,nn_index)= newcentroids_sum(1,nn_index)+r*x;
 newcentroids_sum(2,nn_index)= newcentroids_sum(2,nn_index)+r*y;
 newcentroids_sum(3,nn_index)= newcentroids_sum(3,nn_index)+r;
end

end
function C = cartesian(varargin)
    args = varargin;
    n = nargin;

    [F{1:n}] = ndgrid(args{:});

    for i=n:-1:1
        G(:,i) = F{i}(:);
    end

    C = unique(G , 'rows');
end




function [newcentroids_sum]=zero_sums(newcentroids_sum)

   [x,y]=size(newcentroids_sum);
   for i=1:x
       for  j=1:y
       newcentroids_sum(i,j)=0;
       end
   end
   

end

function [pixel]= clr_image (height,width,pixel)

for i=1:width
    for j=1:height
        pixel(j,i)= 255;
    end
end
end

function [image]=draw_pointswith_magnify (centroids,h,w)
global magnify;
magnify=8;
magnified_h=h*magnify;
magnified_w=w*magnify;
magnifiedsize=[magnified_h,magnified_w];
whiteimage=255 * ones(h,w,'uint8');
imshow(whiteimage);
pixel=zeros(magnified_h,magnified_w);
for y=1:magnified_h
    for x=1:magnified_w
        
        pixel(y,x)=255;
    end
end
clr_image(magnified_h,magnified_w,pixel);
magnified_centroids=zeros(2,w);
for i=1:w
    for j=1:2
       magnified_centroids(j,i)=magnify*centroids(j,i);
    end
end

image=draw_points(magnified_centroids , pixel); 

end
  

function [image]=draw_points(magnified_centroids,pixel)
[~,c]=size(magnified_centroids);
for i=1:c
  pts=floor(magnified_centroids(:,i));
  a=pts(1,1);
  b=pts(2,1);
  if a==0 && b==0
        continue
  end
  pixel(a,b)= 0;
end
image=pixel;
end





function [lanes,P]=connect_dots_with_lines(nodes)

lanes=[];
[r,c]=size(nodes)
fprintf('finding tsp ..');
[P,~]=tspsearchn(nodes);


path=[];
 for i=1:r
     seq=P(i);
     node=nodes(seq,:);
     path=[path;node];
     
 end
 
for i=1:r
     
     lanes=[lanes;path(i,:)];
end
end
     
function [p,L] = tspsearchn(X,m)

[n,dim] = size(X);
if dim == 2 || dim == 3
 D = distmat(X);
else
    mess = 'First argument must be Nx2, Nx3 or symmetric and nonnegative.';
    error('tspsearch:first',mess)
end

if nargin < 2 || isempty(m)
    m = 1;
elseif ~isscalar(m) || m < 1 || m > n || fix(m) < m
    mess = 'M must be an integer in the range 1 to %d.';
    error('tspsearch:second',mess,n)
end

% Starting points for nearest neighbour tours
s = randperm(n);

Lmin = inf;
for k = 1:m
    % Nearest neighbour tour
	p = greedy(s(k),D);
    % Improve tour by 2-opt heuristics
	[p,L] = exchange2(p,D);
    % Keep best tour
	if L < Lmin
        Lmin = L;
        pmin = p;
    end
end

% Output
p = double(pmin);
L = Lmin;
function D = distmat(X)

[n,dim] = size(X);
D=zeros(n);
for j = 1:n
    for k = 1:dim
        v = X(:,k) - X(j,k);
        D(:,j) = D(:,j) + v.*v;
    end
end
D = sqrt(D);

end

function p = greedy(s,D)

n = size(D,1);
p = zeros(1,n,'uint16');
p(1) = s;

for k = 2:n
    D(s,:) = inf;
    [junk,s] = min(D(:,s)); 
    p(k) = s;
end
end

function [p,L] = exchange2(p,D)

n = numel(p);
zmin = -1;

% Iterate until the tour is 2-optimal
while zmin < 0

    zmin = 0;
    i = 0;
    b = p(n);

    % Loop over all edge pairs (ab,cd)
    while i < n-2
        a = b;
        i = i+1;
        b = p(i);
        Dab = D(a,b);
        j = i+1;
        d = p(j);
        while j < n
            c = d;
            j = j+1;
            d = p(j);
            % Tour length diff z
            % Note: a == d will occur and give z = 0
            z = (D(a,c) - D(c,d)) + D(b,d) - Dab;
            % Keep best exchange
            if z < zmin
                zmin = z;
                imin = i;
                jmin = j;
            end
        end
    end

    % Apply exchange
    if zmin < 0
        p(imin:jmin-1) = p(jmin-1:-1:imin);
    end

end

% Tour length
q = double(p);
ind = sub2ind([n,n],q,[q(2:n),q(1)]);
L = sum(D(ind));
end
end



