function [] = plot_time_dgbsv()
    n = [25 50 100 250 400 600 800 999 1500 2500 5000 10000]
    time_row = [13 28 58 131 215 325 435 545 777 1272 2589 5122]
    time_col = [11 23 43 109 185 272 364 465 653 1092 2192 4307]
    
    //title("Temps d`exécution dgbsv")
    plot(n, time_row, "r", n, time_col, "b")
    legend(['row major';'col major'], opt=2)
    xlabel("n")
    ylabel("temps exécution (µs)")
endfunction

function [] = plot_time_dgbmv()
    n = [25 50 100 250 400 600 800 999 1500 2500 5000 10000]
    time_dgbmv = [0 1 2 6 10 16 21 26 38 64 132 259]
    
    //title("Temps d`exécution dgbsv")
    plot(n, time_dgbmv)
    legend(['dgbmv'], opt=2)
    xlabel("n")
    ylabel("temps exécution (µs)")
endfunction

function [] = plot_time_lu()
    n = [25 50 100 250 400 600 800 999 1500 2500 5000]
    time_lu_band = [0 2 5 14 25 35 47 60 90 154 334]
    time_lu_dense = [1 3 6 22 44 63 91 198 420 950 2585]
    
    //title("Temps d`exécution dgbsv")
    plot(n, time_lu_band, "r", n, time_lu_dense, "b")
    legend(['LU band';'LU dense'], opt=2)
    xlabel("n")
    ylabel("temps exécution (µs)")
endfunction
