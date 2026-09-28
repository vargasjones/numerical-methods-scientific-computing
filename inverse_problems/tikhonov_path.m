function X = tikhonov_path(A, b, lambdas)
%TIKHONOV_PATH Compute zero-order Tikhonov solutions for many parameters.
%
% For each lambda > 0, solve
%   min_x ||Ax-b||_2^2 + lambda ||x||_2^2
% using an economy-size SVD of A.

    if any(lambdas <= 0)
        error('All regularization parameters must be positive.');
    end

    [U,S,V] = svd(A,'econ');
    s = diag(S);
    projected_data = U' * b;

    X = zeros(size(A,2), numel(lambdas));
    for j = 1:numel(lambdas)
        filter = s ./ (s.^2 + lambdas(j));
        X(:,j) = V * (filter .* projected_data);
    end
end
