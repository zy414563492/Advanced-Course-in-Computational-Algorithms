nrows = 1000
ncols = 1000
density = 0.005

A = sprand(nrows, ncols, density)
A = A.' + A
A_D = diag(A)
A_D_matrix = diag(A_D)
A_L = tril(A, -1)

x = rand(nrows, 1)
z1 = (A_L + A_D_matrix + A_L') * x

function[M] = MatVecMultiplicationByCRS(A)
    nnzs = nnz(A)
    [ij,val]=spget(A)
    col_ind = ij(:, 2)
    row_ind = ij(:, 1)
    
    row_ptr = []
    cur_ind = 0
    for k = 1:nnzs
        if cur_ind <> row_ind(k) then
            for i = 1:row_ind(k)-cur_ind
                row_ptr = [row_ptr k]
            end
            cur_ind = row_ind(k)
        end
    end
    
    for i = 1:nrows-nnz(row_ptr)+1
        row_ptr($+1) = nnzs+1
    end
    
    for i = 1:nrows
        M(i) = 0.0
        for j = row_ptr(i):row_ptr(i+1)-1
            M(i) = M(i) + val(j) * x(col_ind(j));
        end
    end
endfunction

for i = 1:nrows
    diagVec(i) = x(i) * A_D(i)
end

z2 = MatVecMultiplicationByCRS(A_L) + diagVec + MatVecMultiplicationByCRS(A_L')
result = norm(z1-z2)
