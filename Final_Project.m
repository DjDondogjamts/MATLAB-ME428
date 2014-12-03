clear all;
close all;
clc;

tic         % start clock

h = 0.1;           % R-K step size
x = 0:h:8;         % calculate x values
eps = 0.00001;      % epsilon, given
e1 = 0;             % initialize error
e2 = 0;             % initialize error
ds = 0.001;          % Newton-Raphson step size

y1 = zeros(1,length(x));     % initialize y1
y2 = zeros(1,length(x));     % initialize y2
y3 = zeros(1,length(x));     % initialize y3

y1(1) = 0;                  % given
y2(1) = 0;                  % given
y2(length(x)) = 1;          % given

s = 0.1;                   % initial guess for y3(1)

k = 0;                      % variable to count number of iterations required to satisfy epsilon


for j = 1:100        % If too many iterations, stop and change initial guess
    
    k = k + 1;              % count number of iterations performed
    
    y3(1) = s;              % set initial guess
    
    for i=1:(length(x)-1)
    
        % Solve for Runge-Kutta Constants and calculate again
    
        k1_y1 = y2(i) * h;
        k1_y2 = y3(i) * h;
        k1_y3 = -y1(i) * y3(i) * h;
    
        k2_y1 = (y2(i) + k1_y1 / 2) * h;
        k2_y2 = (y3(i) + k1_y2 / 2) * h;
        k2_y3 = -(y1(i) + k1_y1 / 2) * (y3(i) + k1_y3 / 2) * h;
    
        k3_y1 = (y2(i) + k2_y1 / 2) * h;
        k3_y2 = (y3(i) + k2_y2 / 2) * h;
        k3_y3 = -(y1(i) + k2_y1/2) * (y3(i) + k2_y3 / 2) * h;
    
        k4_y1 = (y2(i) + k3_y1) * h;
        k4_y2 = (y3(i) + k3_y2) * h;
        k4_y3 = -(y1(i) + k3_y1) * (y3(i) + k3_y3) * h;
    
        y1(i+1) = y1(i) + (1/6) * (k1_y1 + 2*k2_y1 + 2*k3_y1 + k4_y1) * h;
        y2(i+1) = y2(i) + (1/6) * (k1_y2 + 2*k2_y2 + 2*k3_y2 + k4_y2) * h;
        y3(i+1) = y3(i) + (1/6) * (k1_y3 + 2*k2_y3 + 2*k3_y3 + k4_y3) * h;
    
    end
    
    e1 = y2(length(x));     % store final value of first y2
    y3(1) = s + ds;     % add Newton-Raphson step to initial guess and recompute

    for i=1:(length(x)-1)   

        % Solve for Runge-Kutta Constants and calculate again
    
        k1_y1 = y2(i) * h;
        k1_y2 = y3(i) * h;
        k1_y3 = -y1(i) * y3(i) * h;
    
        k2_y1 = (y2(i) + k1_y1 / 2) * h;
        k2_y2 = (y3(i) + k1_y2 / 2) * h;
        k2_y3 = -(y1(i) + k1_y1 / 2) * (y3(i) + k1_y3 / 2) * h;
    
        k3_y1 = (y2(i) + k2_y1 / 2) * h;
        k3_y2 = (y3(i) + k2_y2 / 2) * h;
        k3_y3 = -(y1(i) + k2_y1/2) * (y3(i) + k2_y3 / 2) * h;
    
        k4_y1 = (y2(i) + k3_y1) * h;
        k4_y2 = (y3(i) + k3_y2) * h;
        k4_y3 = -(y1(i) + k3_y1) * (y3(i) + k3_y3) * h;
    
        y1(i+1) = y1(i) + (1/6) * (k1_y1 + 2*k2_y1 + 2*k3_y1 + k4_y1) * h;
        y2(i+1) = y2(i) + (1/6) * (k1_y2 + 2*k2_y2 + 2*k3_y2 + k4_y2) * h;
        y3(i+1) = y3(i) + (1/6) * (k1_y3 + 2*k2_y3 + 2*k3_y3 + k4_y3) * h;
    
    end

    e2 = y2(length(x));     % store final value of second y2
    
    s = s - (e1-1)/((e2-e1)/ds);    % recompute the initial guess using Newton-Raphson method
    
    error = abs(e1-1);  
    
    % Plot current values for y1, y2, and y3
    hold on
    plot(x,y1,'-',x,y2,'--',x,y3,'-.')
    hold off
    
    if error < eps      % check if with pre-determined tolerance
        
        fprintf('Number of Iterations = %d \nFinal Error = %4.15f \n\n',k,error);
        
        % add title, line, and labels
        title({'Determining the Flow of Fluid over a Flat Plate Numerically', 'using 4th Order Runge-Kutta and the Shooting Method'})
        xlabel('\eta')
        line('XData', [0 8], 'YData', [1 1], 'LineStyle', '-', 'LineWidth', 2, 'Color','m')
        legend('F','F''','F''''')
        
        break
        
    end
    
    
    
end


if j > 100
    
    fprintf('Too many iterations required, adjust initial guess and try again\n\n');
    
end



