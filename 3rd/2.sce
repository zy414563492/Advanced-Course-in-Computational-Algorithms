clear
funcprot(0)

n = 6
ij = [2, 1; 3, 1; 4, 1; 1, 3; 4, 3; 5, 3; 2, 4; 5, 4; 6, 4; 6, 5; 4, 6; 5, 6]
g = ones(1, size(ij, 1))
G = sparse(ij, g, [n, n])
G = full(G)
c = sum(G, 1)

p = 0.8
d = (1 - p) / n
A = ones(n, n) / n

for i = 1:n
    for j = 1:n
        if c(j) ~= 0
            A(i, j) = p * G(i, j) / c(j) + d
        end
    end
end

v = rand(n, 1)

// Power method
for k = 1:20
    w = A * v
    lambda = v' * w
    v = w / norm(w)
end
