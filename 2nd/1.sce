clear
funcprot(0)

// get results from CRS & CCS format
function[y] = CRS2Vec(val, col_ind, row_ptr, x)
    for i = 1:length(row_ptr)-1
        y(i) = 0.0
        for j = row_ptr(i):row_ptr(i+1)-1
            y(i) = y(i) + val(j) * x(col_ind(j));
        end
    end
endfunction

function[y] = CCS2Vec(val, row_ind, col_ptr, x)
    n = length(col_ptr)-1
    y = zeros(n, 1)
    for j = 1:n
        for i = col_ptr(j):col_ptr(j+1)-1
            y(row_ind(i)) = y(row_ind(i)) + val(i) * x(j);
        end
    end
endfunction


// CG method
function[x,index] = CG(AD, AL, col_ind, row_ptr, b, x, e, k)
    r = b - (CRS2Vec(AL, col_ind, row_ptr, x) + CCS2Vec(AL, col_ind, row_ptr, x) + AD.*x)
    p = r
    for i = 1:k
        q = CRS2Vec(AL, col_ind, row_ptr, p) + CCS2Vec(AL, col_ind, row_ptr, p) + AD.*p
        Alpha = (r'*r)/(p'*q)
        x = x + Alpha * p
        r_new = r - Alpha * q
        index(i) = sqrt(r_new'*r_new)/sqrt(b'*b)
        if index(i) <= e
            break
        end
        Beta = (r_new'*r_new)/(r'*r)
        p = r_new + Beta * p
        r = r_new
    end
endfunction


// CG method with IC(0) preconditioner
function[L, D] = IC0(AD, AL, col_ind, row_ptr)
    n = length(AD)
    nz = length(AL)
    D = AD
    L = zeros(nz, 1)
    for i = 1:n
        w = zeros(i-1, 1)
        for j = row_ptr(i):row_ptr(i+1)-1
            w(col_ind(j)) = AL(j)
            for k=row_ptr(col_ind(j)):row_ptr(col_ind(j)+1)-1
                w(col_ind(j)) = w(col_ind(j)) - L(k) * w(col_ind(k))
            end
            L(j) = w(col_ind(j))/D(col_ind(j))
        end
        for j=row_ptr(i):row_ptr(i+1)-1
            D(i) = D(i) - L(j) * w(col_ind(j))
        end
    end
endfunction

function z = LDLTsolve(L, D, r, col_ind, row_ptr)
    n = length(D)
    z = r
    for i=1:n
        for j=row_ptr(i):row_ptr(i+1)-1
            z(i) = z(i) - L(j) * z(col_ind(j))
        end
    end
    
    for i=i:n
        z(i) = z(i)/D(i)
    end
    
    for i=n:-1:1
        for j=row_ptr(i+1)-1:-1:row_ptr(i)
            z(col_ind(j)) = z(col_ind(j)) - L(i) * z(i)
        end
    end
endfunction

function[x,index] = CGIC0(L, D, AD, AL, col_ind, row_ptr, b, x, e, k)
    r = b - (CRS2Vec(AL, col_ind, row_ptr, x) + CCS2Vec(AL, col_ind, row_ptr,x) + AD.*x)
    z = LDLTsolve(L, D, r, col_ind, row_ptr)
    p = z
    for i = 1:k
        q = CRS2Vec(AL, col_ind, row_ptr, p) + CCS2Vec(AL, col_ind, row_ptr, p) + AD.*p
        a = (r'*z)/(p'*q)
        x = x + a*p
        r_new = r - a*q
        index(i) = sqrt(r_new'*r_new)/sqrt(b'*b)
        if index(i) <= e
            break
        end
        z_new = LDLTsolve(L, D, r_new, col_ind, row_ptr)
        Beta = (r_new'*z_new)/(r'*z)
        p = z_new + Beta*p
        z = z_new
        r = r_new
    end
endfunction


exec('GenLS.sci');

density = 0.005
eTOL = 10^-12;
k = 400 // roop times

for N = 101 //11:10:101
    x = sprand((N-1)**2, 1, density)
    [AD, AL, col_ind, row_ptr, b] = GenLS(N)
    
    [x_CG, index_CG] = CG(AD, AL, col_ind, row_ptr, b, x, eTOL, k)
    [L, D] = IC0(AD, AL, col_ind, row_ptr)
    [x_CGIC0, index_CGIC0] = CGIC0(L, D, AD, AL, col_ind, row_ptr, b, x, eTOL, k)
    
    // plot
    figure(1)
    clf
    xlabel("k")
    ylabel("index: log(y)")
    title("N: " + string(N))
    plot(log10(index_CG), 'r')
    plot(log10(index_CGIC0), 'b')
    legend("index_CG", "index_CGIC0");
    disp("x_CG:", x_CG)
    disp("x_CGIC0:", x_CGIC0)
end

