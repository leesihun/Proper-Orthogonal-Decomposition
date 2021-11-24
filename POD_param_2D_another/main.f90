
    program Main
    use var_POD
    implicit none
    double precision r8_epsilon
    !external r8_f2
    integer i, itime, ntime
    double precision t_start, t_stop
    double precision t, t_out
    double precision work(100+21*Nmodes), yout(Nmodes), dydx(Nmodes)
    double precision h
    
    double precision relerr, abserr
    integer iwork(5), flag
    
    double precision start, finish
    
    
    
    
    
    
    
    allocate(AVGu(nnode), AVGv(nnode), AVGp(nnode), Xpos(nnode), Ypos(nnode))
    allocate(PODu(nnode, Nmodes), PODv(nnode, Nmodes), PODp(nnode, Nmodes))
    allocate(Gvel_u(nnode), Gvel_v(nnode)); Gvel_u=0.d0; Gvel_v=0.d0
    allocate(A(Nmodes), B(Nmodes, Nmodes), C(Nmodes, Nmodes, Nmodes))
    A=0.d0; B=0.d0; C=0.d0
    
    allocate(Ag(Nmodes), Bg(Nmodes, Nmodes)); Ag=0.d0; Bg=0.d0
    
    allocate(D(Nmodes, Nmodes)); D=0.d0
    allocate(E(Nmodes)); E=0.d0
    allocate(F(Nmodes, Nmodes)); F=0.d0
    allocate(G(Nmodes, Nmodes, Nmodes)); G=0.d0
    allocate(KMg(Nmodes)); KMg=0.d0
    
    allocate(HM(Nmodes, Ncntrl)); HM=0.d0
    allocate(JM(Nmodes, Nmodes, Ncntrl)); JM=0.d0
    allocate(KM(Nmodes, Ncntrl)); KM=0.d0
    allocate(LM(Nmodes, Ncntrl, Ncntrl)); LM=0.d0
    allocate(LMg(Nmodes, Ncntrl)); LMg=0.d0
    
    
    
    call cpu_time(start)
    
    !call Re_POD_new2(nnode, nelem)
    !stop
    !call Get_POD_new(nnode, nelem)
    !call POD_cfd_cntrl (1)
    !call Get_POD_new2(nnode, nelem)
    !call Get_POD_new4(nnode, nelem)
    !call Get_POD_new6(nnode, nelem)
    
    
    
    
    !!! USE THIS to get reconstructed flow
    !call Get_POD_new8_overset2(nnode, nelem)
    
    !ONLY POD MODES AND COEFFS.
    call Get_POD_new8_overset4(nnode, nelem)
    
    
    
    
    
    
    
    
    
    
    !call Get_POD_new8(nnode, nelem)
    !call POD_cfd_ale (1)
    !call POD_cfd_cntrl (1)
    
    
    
    
        ! put code to test here
    call cpu_time(finish)
    print*, 'Elapsed Time = ',finish-start
    
    
    pause;
    pause;
    stop
    
    
    !! Coefficients, normally on
    !call POD_cfd_cntrl_ale (1)
    
    stop
    
    pause;
    pause;
    
    end program Main
    
    