function problem4

%% Konfiguration
h = 0.1;
dt = 0.001;


%% Startwerte

data = double(imread('problem4-2.tiff'));

[m,n] = size(data);

[x,y] = meshgrid(0:h:(m-1)*h,0:h:(n-1)*h);

u = ones(m,n);

index_x = 1:n;
index_y = 1:m;

center_x = round(length(index_x)/2);
center_y = round(length(index_y)/2);
blockSize = 15;
block_x = center_x-blockSize:center_x+blockSize;
block_y = center_y-blockSize:center_y+blockSize;

% Temperaturleitungskoeffizient für das Gitter:
c = 1.25 * data/255;

% Simulations-figure
figure('Position', [0, 0, 1000, 600]);
handle = surf(x,y,u);
axis manual;
%axis([x_min-border, x_max+border, y_min-border, y_max+border, -1, 1]);
axis([0, h*n, 0, h*m, -1, 5]);

drawnow

%% Simulation
t = 0;

count = 0;

while 1
    u = dt*4*c.*del2(u,h)+u;
    u(block_x, block_y) = u(block_x, block_y) + 10*dt;
    
    % Neumannrand:
    u(:,1) = u(:,2);
    u(1,:) = u(2,:);
    u(m,:) = u(m-1,:);
    u(:,n) = u(:,n-1);

%     % Dirichletrand:
%     u(:,1) = 0;
%     u(1,:) = 0;
%     u(mx,:) = 0;
%     u(:,nx) = 0;

    t = t + dt;
    count = count+1;

    if(mod(count, 1000) == 0) 
        set(handle, 'ZData', log(u));
        drawnow
        
        fprintf('t = %f\n', t);
    end
end


end

