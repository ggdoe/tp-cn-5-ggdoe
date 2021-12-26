exec lu-tridiag.sci;

function [x] = it_jacobi(x, A, b)
    x = x + 1/2 * (b-A*x)
endfunction

function [x] = it_gauss_seidel_vec(x, A, b)
    // D-E = tril(A) ; F = tril(A)-A
    x = tril(A)\((tril(A)-A)*x + b)
endfunction

function [err] = jacobi_poisson1D(n, nbr_it)
    A = make_tridiag(n, -1, 2, -1)
    x = zeros(n, 1)
    b = get_RHS(n, -5, 5)
    ex = get_analytic_sol(n, -5, 5)
    err = zeros(nbr_it, 1)
    err(1) = norm(x-ex)/norm(ex)
    
    for i=2:nbr_it
        x = it_jacobi(x, A, b)
        err(i) = norm(x-ex)/norm(ex)
    end
endfunction

function [err] = seidel_poisson1D(n, nbr_it)
    A = make_tridiag(n, -1, 2, -1)
    x = zeros(n, 1)
    b = get_RHS(n, -5, 5)
    ex = get_analytic_sol(n, -5, 5)
    err = zeros(nbr_it, 1)
    err(1) = norm(x-ex)/norm(ex)
    
    for i=2:nbr_it
        x = it_gauss_seidel_vec(x, A, b)
        err(i) = norm(x-ex)/norm(ex)
    end
endfunction

function [] = plot_jac_seid(n, nbr_it)
    f = scf(1);
    f.color_map = rainbowcolormap(32);
    
    err = jacobi_poisson1D(n, nbr_it)
    plot2d(1:nbr_it, err, style=1, leg="Jacobi")
    
    err = seidel_poisson1D(n, nbr_it)
    plot2d(1:nbr_it, err, style=15, leg="Gauss-Seidel")
    
    //title("Temps d`exécution dgbsv")
    legend(["Jacobi" "Gauss-Seidel"], opt=3)
    xlabel("nombre itération")
    ylabel("erreur relative")
endfunction

function [] = time_jac_seid()
    n = [10 20 30 40 50 60 70 80 90 100 120 140 160 180 200]
    size_n = size(n)(2)
    f = scf(1); f.color_map = rainbowcolormap(32);
    time_j = zeros(size(n)(2),1)
    time_gs = zeros(size(n)(2),1)
    nbr_rep = 10 ; eps = 0.05
    
  for i = 1:size_n
        A = make_tridiag(n(i), -1, 2, -1)
        b = get_RHS(n(i), -5, 5)
        ex = get_analytic_sol(n(i), -5, 5)
        
        for j = 1:nbr_rep
            x = zeros(n(i), 1)
            err = 1
            tic()
            while err > eps
                x = it_jacobi(x, A, b)
                err = norm(x-ex)/norm(ex)
            end
            time_j(i) = time_j(i) + toc()
            
            x = zeros(n(i), 1)
            err = 1
            tic()
            while err > eps
                x = it_gauss_seidel_vec(x, A, b)
                err = norm(x-ex)/norm(ex)
            end
        end
        time_gs(i) =  time_gs(i) + toc()
    end
    time_j = time_j / nbr_rep
    time_gs = time_gs / nbr_rep
    plot2d(n, time_j, style=1, leg="Jacobi")
    plot2d(n, time_gs, style=15, leg="Gauss-Seidel")
    
    //title("Temps d`exécution dgbsv")
    legend(["Jacobi" "Gauss-Seidel"], opt=1)
    xlabel("nombre itération")
    ylabel("Temps d`éxécution (s)")
endfunction



function [ex] = get_analytic_sol(n, T0, T1)
    ex  = zeros(n, 1)
    h = 1/(n+1)
    for i = 1:n
        ex(i) = T0 + i*h * (T1 - T0)
    end
endfunction

function [b] = get_RHS(n, T0, T1)
    b = zeros(n, 1)
    b(1) = T0
    b(n) = T1
endfunction
