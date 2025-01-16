function T = chebyshev_polynomials(x, k_nodes)
    n = length(k_nodes);
    T = zeros(n, 1);
    T(1) = 1;
    T(2) = x;

    for i = 3:n
        T(i) = 2 * x * T(i - 1) - T(i - 2);
    end
end
