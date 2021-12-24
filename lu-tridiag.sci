function [A] = make_tridiag(n, a, b, c)
    A = a*diag(ones(n-1,1),1) + b*eye(n,n) + c*diag(ones(n-1,1),-1)
endfunction

function [A] = make_tridiag_vec(a, b, c)
    A = diag(a,1) + diag(b) + diag(c,-1)
endfunction
