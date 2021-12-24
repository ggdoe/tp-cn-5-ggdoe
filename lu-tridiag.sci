function [A] = make_tridiag(n, a, b, c)
    A = b*eye(n,n) + a*diag(ones(n-1,1),1) + c*diag(ones(n-1,1),-1)
endfunction
