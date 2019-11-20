clear
nrows = 1000
ncols = 1000
density = 0.005
A =sprand(nrows, ncols, density)
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


//disp(val)
//disp(col_ind)
//disp(row_ptr)

