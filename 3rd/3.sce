clear
funcprot(0)

//A_ori = [165 158 153 174 171 157 177 163 164 172; 67 56 48 68 62 49 79 56 58 70]'
A_ori = strtod(read_csv('./femail_middel.csv')(2:37, 2:3))
A = A_ori
m1 = mean(A(:, 1))
m2 = mean(A(:, 2))
A(:, 1) = A(:, 1) - m1
A(:, 2) = A(:, 2) - m2
C = A' * A

v = rand(2, 1)

// Power method
for k = 1:20
    w = C * v
    lambda = v' * w
    v = w / norm(w)
end

v = v / norm(v)
u = A * v
u = u / norm(u)
lambda = v' * (C * v) / (v' * v)
s = sqrt(lambda)
P = s * u * v'

P(:, 1) = P(:, 1) + m1
P(:, 2) = P(:, 2) + m2

disp(P(:, 1))

// plot
figure(1)
clf
xlabel("Height")
ylabel("Weight")
scatter(A_ori(:, 1), A_ori(:, 2), 'red')
scatter(P(:, 1), P(:, 2), 'blue')
plot(P(:, 1), P(:, 2), 'b')
legend("A", "P");



