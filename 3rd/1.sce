clear
funcprot(0)

A = [0.05 0.85 0.05 0.25; 0.45 0.05 0.45 0.25; 0.45 0.05 0.05 0.25; 0.05 0.05 0.45 0.25]
v = rand(4, 1)
res = [0.62427; 0.57716; 0.41225; 0.32745]

// Power method
for k = 1:20
    w = A * v
    lambda = v' * w
    v = w / norm(w)
end

check = norm(v - res)
