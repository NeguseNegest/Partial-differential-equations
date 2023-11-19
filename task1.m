
clc




% Generate 1000 points in [0,1]x[0,1]
M = 1000;
xpoints = rand(M,1);
ypoints = rand(M,1);
k = 2*pi;
N = 7;

% Initialize v1 and v2
v1 = ypoints;  % v1 = y
v2 = 1 - xpoints;  % v2 = 1 - x
%f_old=zeros(1,length(v1));
% Initialize matrix A
A = zeros(M, (2*N + 1)^2);
f_old = zeros(M, 1); 
% Populate the matrix A
for i = 1:length(xpoints)
    xp = xpoints(i);
    yp = ypoints(i);
    for n1 = -N:N
        for n2 = -N:N
            fourier_term = exp(1i * k * (n1 * xp + n2 * yp));
            col = (n1 + N) * (2*N + 1) + (n2 + N) + 1; % Find the column corresponding to n1 and n2
            A(i, col) = fourier_term;
        end
    end
end

% TODO: 
% Populate matrix A - done
% Populate v1 and v2 - done
% Need to solve for v_new: A*A^T*v_new = A^T*v_old
% Need v_new
v_hat_1 = (A' * A) \ (A' * v1);  
v_hat_2=(A' * A) \ (A' * v2);
% Obtain the real part of v_hat_1
v_bar_1 = real(v_hat_1);
v_bar_2=real(v_hat_2);

% --------------------Populate the f vector ------------


for i = 1:length(xpoints)
    if sqrt((xpoints(i) - 0.5)^2 + (ypoints(i) - 0.5)^2) <= 0.1
        f_old(i) = 1;
    end
end
%---------------------end----------------------------------

%------------------Find f_new-------------------
f_im=(A' * A) \ (A' * f_old);
f_new=real(f_im);% Get the real part of the vector 

%-----------------end---------------------------


