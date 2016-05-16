function drawResult()

u_h = textread('output.txt');
n = sqrt(length(u_h));
u_h = reshape(u_h,n,n);
[X,Y] = meshgrid(1:n,1:n);
surf(X,Y,u_h);

end