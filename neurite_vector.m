function [U, V, theta] = neurite_vector(image, sigma)
%%
% INPUT:
% sigma: size of gaussian kernel
% image: image to process

% 0 - 1 double
image = double(image);
image = image ./ max(image(:));

gaussian_img = imgaussfilt(image, sigma, 'padding', 'replicate');

% second-order derivatives
[gx, gy] = gradient(gaussian_img);
[gxx, gxy] = gradient(gx);
[gyx, gyy] = gradient(gy);

[m, n] = size(image);
E = zeros(m, n);
U = E;
V = E;

for i = 1:m
   for j = 1:n 
        Hessian = [gxx(i,j) gxy(i,j); gyx(i,j) gyy(i,j)];
        [eigen_vector, eigen_value] = eig(Hessian, 'vector');
        [~, index] = max(abs(eigen_value));
        lambda_x = eigen_value(index);
        E(i, j) = (lambda_x < 0) * lambda_x;
        % eigenvector
        U(i, j) = eigen_vector(1, 3 - index);
        V(i, j) = eigen_vector(2, 3 - index);
   end
end

% eigenvector .* eigenvalue
U = U .* E;
V = V .* E;

theta = atan2d(V, U);