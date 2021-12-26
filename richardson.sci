exec lu-tridiag.sci;

function [x] = it_richardson(x, alpha, A, b)
    x = x + alpha * (b-A*x)
endfunction

function [err] = richardson_poisson1D(n, nbr_it, alpha)
    A = make_tridiag(n, -1, 2, -1)
    x = zeros(n, 1)
    b = get_RHS(n, -5, 5)
    ex = get_analytic_sol(n, -5, 5)
    err = zeros(nbr_it, 1)
    err(1) = norm(x-ex)/norm(ex)
    
    for i=2:nbr_it
        x = it_richardson(x, alpha, A, b)
        err(i) = norm(x-ex)/norm(ex)
    end
endfunction

function [] = plot_richardson(n, nbr_it)
    alpha = [0.02 0.05 0.1 0.15 0.25 0.35 0.45 0.50 0.5+2/nbr_it]
    f = scf();
    f.color_map = rainbowcolormap(32);
    
    for i = 1:size(alpha)(2)
        err = richardson_poisson1D(n, nbr_it, alpha(i))
        plot2d(1:nbr_it, err, style=3*i, leg=string(alpha(i)))
    end
    //title("Temps d`exécution dgbsv")
    legend("alpha = " + string(alpha), opt=3)
    xlabel("nombre itération")
    ylabel("erreur relative")
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
