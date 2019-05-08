function [z_w varargout] = myWhiten(z)


% Get dimensions
[m n] = size(z);

% Compute sample covariance
R = zeros(m);
mean_z = mean(z,2);
for i = 1:m
    for j = 1:i
        num = 0;
        for k = 1:n
           num = num + (z(i,k) - mean_z(i)) * (z(j,k) - mean_z(j));
        end
        R(i,j) = num / (n-1);
        R(j,i) = R(i,j);
    end
end

% Whiten z
[U D ~] = svd(R,'econ');
DIM  = size(D,1);
T = zeros(m);
for i = 1:DIM
    T = T + 1/sqrt(D(i,i)) * U(:,i) * U(:,i)';
end
z_w = T * z;

if (nargout == 2)
    varargout{1} = T;
end
