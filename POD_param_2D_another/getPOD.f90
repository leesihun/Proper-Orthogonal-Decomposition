    subroutine Get_POD_new8_overset4(m_dim, nel)
    use var_POD
    implicit none
    integer i, i1, j, k, istart, n_snap,n_snap_t, nsnp
    integer m_dim, nel, dt_val,nwall
    double precision, parameter :: dt = 1.d0
    double precision, parameter :: t_start1=1, t_start2=1001, t_start3=300  ! org : t_start2 = 4500
    character(len=255) time, root(Ncntrl+2), char1, param1, dir1,param2
    double precision, parameter :: Uref=1.0d0   ! reference speed
    double precision bi,bj,bk, ci,cj,ck, area, t_start
    double precision xel(3), yel(3), gam_val, pi
    integer parameters
    
    double precision  t_read, tmpppp
    double precision, allocatable :: u_avg(:), v_avg(:), eval(:)
    double precision, allocatable :: x(:), y(:) !, conv(:,:)
    double precision, allocatable :: u(:,:), v(:,:)
    double precision, allocatable :: u1(:,:), v1(:,:)
    double precision, allocatable :: u_vec(:,:), u_f(:), tmp(:,:), Fcoeff(:,:), temp_f(:)
    double precision, allocatable :: re_u(:,:), re_v(:,:)
    double precision, allocatable :: gam(:), ConFM(:,:), u_avg1(:,:), v_avg1(:,:),gam2(:)
    double precision finish1, start1, finish2, finish3, finish4, finish5
    
    double precision cht_l, cht_u, nd_freq, phy_freq,dt_real
    double precision nd_freq1, phy_freq1, dt_real1, nd_freq2(2)
    character(len=255) root1
    double precision, allocatable :: u2(:,:), v2(:,:), Fcoeff2(:,:)
    double precision, allocatable :: wall_v(:,:),wall_v1(:,:),wall_v2(:,:)
    integer,allocatable :: wallid(:)
    
    
        integer :: Nparam=7
        integer :: param(7)
        
        param(1) = 100
        param(2) = 150
        param(3) = 200
        param(4) = 250
        param(5) = 300
        param(6) = 350
        param(7) = 400
        
        
    dt_val=1;
    time="a"
    
    
    call cpu_time(start1)
    
    !!!!!!nd_freq=0.0
    dt_real=0.01
    
    !!!!!!nd_freq1=0.0
    dt_real1=0.01
    
    !dt_real1=1.396263402401660e-2
    PI=4.D0*DATAN(1.D0)
    !write(*,*) dt_real1,dt_real,dt_real*300,dt_real1*3*300
    !pause
    !pause
    !dt_real=0.2d0
    
    !!!!!!dt_real=dt_real*dble(dt_val)
    !!!!!!dt_real1=dt_real1*dble(dt_val)
    
    nd_freq2(1)=0; nd_freq2(2)=0
    
    
    
    t_start=t_start1
    !!!!!!n_snap=n_snapshot*2
    n_snap = n_snapshot ! Newly added... 20/08/27
    n_snap_t=n_snap*Nparam
    
    !!!!!!allocate(wall_v(m_dim,2)); wall_v=0.d0
    !!!!!!allocate(wall_v1(m_dim,2)); wall_v1=0.d0
    !!!!!!allocate(wall_v2(m_dim*(Ncntrl+1),2) ); wall_v2=0.d0
    
    !!!!!!wall_v = 0.d0
    !!!!!!wall_v1 = 0.d0
    !!!!!!wall_v2=0.D0
    !!!!!!wall_v2=0.D0
    !!!!!!allocate(gam(n_snap_t), gam2(n_snap_t))
    !!!!!!gam=0.d0
    !!!!!!do i=1, n_snapshot
        !!!!!!gam(i)=0.d0 !*0.d0
        !!!!!!gam(n_snapshot+i)=0.d0
    !!!!!!enddo
    !!!!!!gam2=gam    !0.d0
    
    
    allocate(gam(n_snap_t), gam2(n_snap_t))
    gam=0.d0
    do i=1, n_snap_t
        gam(i)=0.d0 !*0.d0
        gam(n_snap_t)=0.d0
    enddo
    gam2=gam    !0.d0
    !cht_u=0.1d0
    !cht_l=0.09d0
    !nd_freq=0.d0
    !dt_real=1.d0/0.1957/50.d0
    
    root(1)=adjustl('NACA_plg/')  ! root directory
    
    write(param2, *) param(1)
    param1=adjustl(param2)
    
    dir1 = adjustl('/tecplot_') 
    
    
    !root(1)=adjustl('NACA_plg/1.5/tecplot_')  ! root directory
    
    char1=root(1)
    allocate(weight(m_dim,1)); weight=0.d0
    allocate(x(m_dim), y(m_dim),conv(nel,3))
    write(time,*) floor(t_start)
    time=adjustl(time)
    
    open(1, file=char1(1:len_trim(char1))//param1(1:len_trim(param1))//dir1(1:len_trim(dir1))//time(1:len_trim(time))//'.dat') ! Heaving airfoil
    
    read(1,*)
    read(1,*)
    
    do j=1, m_dim
        read(1,*) x(j), y(j), tmpppp, tmpppp, tmpppp
    enddo
    
    
    do j=1, nel
        read(1,*) conv(j,1), conv(j,2), conv(j,3)
    enddo
    close(1)
    
    !y=y-0.175d0
    
    do i1=1, nel
        do j=1, 3
            xel(j)=x(conv(i1,j))
            yel(j)=y(conv(i1,j))
        enddo
        bi=yel(2)-yel(3)
        bj=yel(3)-yel(1)
        bk=yel(1)-yel(2)
        
        ci=xel(3)-xel(2)
        cj=xel(1)-xel(3)
        ck=xel(2)-xel(1)
        
        area=abs(0.5d0*(-bk*cj+ck*bj));
        do j=1, 3
            weight(conv(i1,j),1)=weight(conv(i1,j),1)+1.d0/3.d0*area
        enddo
    enddo
    
    
    ! Weight done
    
    
    allocate(ia(m_dim,1)); ia=0
    do i=1, nel
        do j=1, 3
            k=conv(i,j);
            ia(k,1)=ia(k,1)+1
        enddo
    enddo
    do i=1, m_dim
        weight(i,1)=weight(i,1)/dble(ia(i,1))
    enddo
    
    allocate(u(m_dim,n_snap_t),v(m_dim,n_snap_t)); u=0.d0; v=0.d0
    allocate(ConFM(m_dim,2), u_avg1(m_dim,Ncntrl+2), v_avg1(m_dim,Ncntrl+2))
    allocate(u_avg(m_dim), v_avg(m_dim)); u_avg=0.d0; v_avg=0.d0
    u_avg1=0.d0; v_avg1=0.d0
    
    !t_start(1) = t_start(1)-dt_val(1);t_start(2) = t_start(2)-dt_val(1);t_start(3) = t_start(2)-dt_val(1)
    
    k=1;
    
    
    
    i1=0; nsnp=0
    
    do parameters = 1, Nparam
        write(*,*) 'Reading parameter' ,parameters
        do i=1, n_snap
            
            i1=i1+1
            gam(i1)=gam_val
        
            t_read=t_start+dt_val*dble(i)
            write(time,*) floor(t_read)
            time=adjustl(time)
            char1=root(1)
            write(*,*) char1(1:len_trim(char1)),time(1:len_trim(time))
            
            write(param2, *) param(parameters)
            param1=adjustl(param2)
            
            open(1, file=char1(1:len_trim(char1))//param1(1:len_trim(param1))//dir1(1:len_trim(dir1))//time(1:len_trim(time))//'.dat') ! Heaving airfoil
            !pause
            read(1,*)
            read(1,*)
        
            do j=1, m_dim
                read(1,*) tmpppp, tmpppp, u(j,i1), v(j,i1), tmpppp          ! read u, v velocity
                u(j,i1) = u(j,i1) 
                v(j,i1) = v(j,i1) 
                u_avg1(j,k)=u_avg1(j,k)+u(j,i1)
                v_avg1(j,k)=v_avg1( j,k)+v(j,i1)
            enddo
            do j=1, nel
                read(1,*) conv(j,1), conv(j,2), conv(j,3)
            enddo
            close(1)
        
        enddo
    enddo
    
    
        
        u_avg1(:,k)=u_avg1(:,k)/dble(n_snap_t)
        v_avg1(:,k)=v_avg1(:,k)/dble(n_snap_t)
        
        u_avg=u_avg+u_avg1(:,k)
        v_avg=v_avg+v_avg1(:,k)
        nsnp=nsnp+n_snap
    
    v_avg=v_avg1(:,1)
    u_avg=u_avg1(:,1)
    
    
    
    print*, 'Reading files complete'
    
    call cpu_time(finish1)
    print*, 'Elapsed Time for reading input files = ',finish1-start1
    
    ConFM(:,1)=0.d0
    ConFM(:,2)=0.d0
    
    open(1,file='avg1.dat')
    write(1,49)
    write(1,47) 1, m_dim, nel
    do j=1, m_dim
        write(1,48) x(j), y(j), &
            u_avg1(j,1), v_avg1(j,1), 0.0
    enddo
    do j=1, nel
        write(1,*) conv(j,1), conv(j,2), conv(j,3)
    enddo
    
    print*, 'avg1.dat written'
    
    allocate(u1(m_dim,n_snap_t),v1(m_dim,n_snap_t)); u1=0.d0; v1=0.d0
    
    do i=1, n_snap_t
        u1(:,n_snap_t*0+i)=u(:,n_snap_t*0+i)-u_avg !-gam(n_snapshot*0+i)*ConFM(:,1)!-gam2(n_snapshot*0+i)*wall_v(:,1)   ! perturbation
        v1(:,n_snap_t*0+i)=v(:,n_snap_t*0+i)-v_avg! -gam(n_snapshot*0+i)*ConFM(:,2)!-gam2(n_snapshot*0+i)*wall_v(:,2)   ! perturbation

    enddo
    
    
    
    deallocate(u, v)
    allocate(u(m_dim, n_snap_t), v(m_dim, n_snap_t))
    allocate(tmp(n_snap_t,n_snap_t))
    allocate(u_vec(n_snap_t, n_snap_t), u_f(n_snap_t))
    
    tmp=0.d0
    do i = 1,n_snap_t
        do j = 1,n_snap_t
            do k=1, m_dim
                tmp(i,j)=tmp(i,j)+weight(k,1)*(u1(k,i)*u1(k,j)+v1(k,i)*v1(k,j))
            enddo
        enddo
    enddo
    
    do i=1, n_snap_t
        do j=1, n_snap_t
            tmp(j,i)=tmp(i,j)
        enddo
    enddo
    tmp=tmp/dble(n_snap_t)
    
    print*, 'start eigen decomposition'
    
    
    
    ! u_f is negative : 순서 바꾸기?
    
    
    
    
    
    call eig_mat(n_snap_t,n_snap_t, tmp, u_vec, u_f)
    deallocate(tmp)
    
    !u_f = -(u_f);!!!!!!!!!!! not original
    
    call cpu_time(finish2)
    print*, 'Elapsed Time for eig_decomposition = ',finish2-finish1
        
    
    
    open(1, file='PODmode/u_f.dat')
    
    print*, 'write POD mode'
    
    
    
    
    
    
    !!!!!!! BELOW ARE ADDED
    !allocate (temp_f(n_snap_t));temp_f = 0.d0
    !do i=1, n_snap_t
    !    !write(1, *) u_f(i)
    !    temp_f(i) = u_f(n_snap_t-(i-1))
    !enddo
    !u_f = temp_f
    !deallocate(temp_f)
    
    
    
    
    
    
    
    
    
    
    
    do i=1, n_snap_t
        write(1,*) u_f(n_snap_t-(i-1))             !%%%%%%%%%%%%% Original code
    enddo
    
    
    close(1)
    
    allocate(tmp(n_snap_t, n_snap_t)); tmp=0.d0
    do i=1, n_snap_t
        tmp(:,i)=u_vec(:,n_snap_t-(i-1))             !%%%%%%%%%%%%% Original code
        
        !tmp(:,i) = u_vec(:,i)
    enddo
    
    
    u_vec=tmp
    deallocate(tmp)
    
    call MM(m_dim, n_snap_t, n_snap_t, n_snap_t, u1, u_vec, u)
    call MM(m_dim, n_snap_t, n_snap_t, n_snap_t, v1, u_vec, v)
    
    !print*, u_f
    !pause;
    
    do i=1, n_snap_t
    !    u(:,i)=u(:,i) /dsqrt(u_f(i)*dble(n_snap_t))
    !    v(:,i)=v(:,i) /dsqrt(u_f(i)*dble(n_snap_t))
        
        u(:,i)=u(:,i) /dsqrt(u_f(n_snap_t-(i-1))*dble(n_snap_t))   !original
        v(:,i)=v(:,i) /dsqrt(u_f(n_snap_t-(i-1))*dble(n_snap_t))   !original
    enddo
    
    do i1=1, Nmodes
        t_read=dble(i1)
        write(time,*) floor(t_read)
        time=adjustl(time)
        char1=root(1)
        write(*,*) 'write PODmode',i1
            
        open(1, file='PODmode/PODmode_'//time(1:len_trim(time))//'.dat') !Heaving/Pitching airfoil
        !pause
        write(1,49)
        write(1,47) floor(t_read), m_dim, nel
        
        do j=1, m_dim
            write(1,48) x(j), y(j), u(j,i1), v(j,i1), 0.d0
        enddo
        do j=1, nel
            write(1,*) conv(j,1), conv(j,2), conv(j,3)
        enddo
        close(1)
    enddo
    
    allocate(Fcoeff(n_snap_t,Nmodes)); Fcoeff=0.d0
    do i=1, Nmodes
        do j=1, n_snap_t
            do k=1, m_dim
                !Fcoeff(j,i)=u_vec(j,i) *dsqrt((u_f(n_snap_t-(i-1)))*dble(n_snap_t))
                Fcoeff(j,i)=Fcoeff(j,i)+weight(k,1)*(u1(k,j)*u(k,i)+v1(k,j)*v(k,i))
            enddo
        enddo
    enddo
    
    open(1,file='Fcoeff_velocity0')
    do i=1, n_snap_t
        write(1,449)(Fcoeff(i,j), j=1,Nmodes)
    enddo
    close(1)
    
   call cpu_time(finish3)
    print*, 'Elapsed Time for POD+coeff. write = ',finish3-finish2 
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !allocate(u2(m_dim,n_snapshot), v2(m_dim,n_snapshot))
    !i1=0
     !   do i=1, n_snapshot !*2
            
         !   i1=i1+1
            !gam(i1)=gam_val(k)
        
        !    t_read=t_start1+dt_val*dble(i)
       !     write(time,*) floor(t_read)
      !      time=adjustl(time)
     !       char1=root(1)
    !        write(*,*) char1(1:len_trim(char1)),time(1:len_trim(time))
            
    !        open(1, file=char1(1:len_trim(char1))//time(1:len_trim(time))//'.dat') !Heaving/Pitching airfoil
            !pause
    !        read(1,*)
    !        read(1,*)
        
    !        do j=1, m_dim
    !            read(1,*) tmpppp, tmpppp, u2(j,i1), v2(j,i1), tmpppp
    !        enddo
    !        do j=1, nel
    !            read(1,*) conv(j,1), conv(j,2), conv(j,3)
    !        enddo
    !        close(1)
        
    !    enddo
    !do i=1, n_snapshot
    !    u2(:,i)=u2(:,i)-u_avg-gam(n_snapshot*0+i)*ConFM(:,1)-gam2(n_snapshot*0+i)*wall_v(:,1)
    !    v2(:,i)=v2(:,i)-v_avg-gam(n_snapshot*0+i)*ConFM(:,2)-gam2(n_snapshot*0+i)*wall_v(:,2)
    !    !write(*,*) u1(1,i),v1(1,i)
    !    !pause
    !enddo   
    !allocate(Fcoeff2(n_snapshot,Nmodes)); Fcoeff2=0.d0
    !do i=1, Nmodes
    !    do j=1, n_snapshot
    !        do k=1, m_dim
    !           !Fcoeff(j,i)=u_vec(j,i) *dsqrt((u_f(n_snap_t-(i-1)))*dble(n_snap_t))
    !            Fcoeff2(j,i)=Fcoeff2(j,i)+weight(k,1)*(u2(k,j)*u(k,i)+v2(k,j)*v(k,i))
    !        enddo
    !    enddo
    !enddo
    !open(1,file='Fcoeff_velocity2_0')
    !do i=1, n_snapshot !*2
    !    write(1,449)(Fcoeff2(i,j), j=1,Nmodes)
    !enddo
    !close(1)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    
    
    
    
    
    
    
    
    allocate(tmp(Nmodes,n_snap_t)); tmp=0.d0
    do i=1, n_snap_t
        do j=1, Nmodes
            tmp(j,i)=Fcoeff(i,j)
        enddo
    enddo
    allocate(re_u(m_dim,n_snap_t), re_v(m_dim,n_snap_t))
    call MM(m_dim, Nmodes, Nmodes, n_snap_t, u, tmp, re_u)
    call MM(m_dim, Nmodes, Nmodes, n_snap_t, v, tmp, re_v)
    deallocate(tmp)
    do i=1, n_snap_t
        re_u(:,i)=re_u(:,i)+u_avg
        re_v(:,i)=re_v(:,i)+v_avg
    enddo
    
    call cpu_time(finish4)
    print*, 'reconstructed vel compute = ',finish4-finish3 
    
    i1=0
    do k=1, nparam
        do i=1, n_snap
            i1=i1+1
            t_read=t_start1+dt*dble(i)
            write(time,*) floor(t_read)
            time=adjustl(time)
            
            
            write(param2, *) param(k)
            param1=adjustl(param2)
            
            !open(1, file=char1(1:len_trim(char1))//param1(1:len_trim(param1))//dir1(1:len_trim(dir1))//time(1:len_trim(time))//'.dat') ! Heaving airfoil
        
            open(1,file='re_construction0/'//param1(1:len_trim(param1))//'/re_'//time(1:len_trim(time))//'.dat')
            write(1,49)
            write(1,47) floor(t_read), m_dim, nel
            do j=1, m_dim
                write(1,48) x(j), y(j),  & !+0.175d0*dcos(2.d0*pi*nd_freq*dt_real*dble(i)), &
                    re_u(j,i1)*Uref, re_v(j,i1)*Uref, 0.d0
            enddo
            
            do j=1, nel
                write(1,*) conv(j,1), conv(j,2), conv(j,3)
            enddo
            close(1)
        enddo
    enddo
    
    call cpu_time(finish5)
    print*, 'write reconstructed vel = ',finish5-finish4
    
    
    !!!!!!deallocate(re_u,re_v,Fcoeff2,u2,v2)
    !!!!!!
    !!!!!!allocate(u2(m_dim,n_snapshot), v2(m_dim,n_snapshot))
    !!!!!!i1=0
    !!!!!!    do i=1, n_snapshot !*2
    !!!!!!        
    !!!!!!        i1=i1+1
    !!!!!!        !gam(i1)=gam_val(k)
    !!!!!!    
    !!!!!!        t_read=t_start2+dt_val(2)*dble(i)
    !!!!!!        write(time,*) floor(t_read)
    !!!!!!        time=adjustl(time)
    !!!!!!        char1=root(3)
    !!!!!!        write(*,*) char1(1:len_trim(char1)),time(1:len_trim(time))
    !!!!!!        
    !!!!!!        open(1, file=char1(1:len_trim(char1))//time(1:len_trim(time))//'.dat') !Heaving/Pitching airfoil
    !!!!!!        !pause
    !!!!!!        read(1,*)
    !!!!!!        read(1,*)
    !!!!!!    
    !!!!!!        do j=1, m_dim
    !!!!!!            read(1,*) tmpppp, tmpppp, u2(j,i1), v2(j,i1), tmpppp
    !!!!!!        enddo
    !!!!!!        do j=1, nel
    !!!!!!            read(1,*) conv(j,1), conv(j,2), conv(j,3)
    !!!!!!        enddo
    !!!!!!        close(1)
    !!!!!!    
    !!!!!!    enddo
    !!!!!!do i=1, n_snapshot
    !!!!!!    u2(:,i)=u2(:,i)-u_avg-gam(n_snapshot*1+i)*ConFM(:,1)-gam2(n_snapshot*1+i)*wall_v1(:,1)
    !!!!!!    v2(:,i)=v2(:,i)-v_avg-gam(n_snapshot*1+i)*ConFM(:,2)-gam2(n_snapshot*1+i)*wall_v1(:,2)
    !!!!!!    !write(*,*) u1(1,i),v1(1,i)
    !!!!!!    !pause
    !!!!!!enddo 
    !!!!!!
    !!!!!!allocate(Fcoeff2(n_snapshot,Nmodes)); Fcoeff2=0.d0
    !!!!!!do i=1, Nmodes
    !!!!!!    do j=1, n_snapshot
    !!!!!!        do k=1, m_dim
    !!!!!!            !Fcoeff(j,i)=u_vec(j,i) *dsqrt((u_f(n_snap_t-(i-1)))*dble(n_snap_t))
    !!!!!!            Fcoeff2(j,i)=Fcoeff2(j,i)+weight(k,1)*(u2(k,j)*u(k,i)+v2(k,j)*v(k,i))
    !!!!!!        enddo
    !!!!!!    enddo
    !!!!!!enddo
    !!!!!!open(1,file='Fcoeff_velocity2_1')
    !!!!!!do i=1, n_snapshot !*2
    !!!!!!    write(1,449)(Fcoeff2(i,j), j=1,Nmodes)
    !!!!!!enddo
    !!!!!!close(1)
    !!!!!!
    !!!!!!allocate(tmp(Nmodes,n_snapshot)); tmp=0.d0
    !!!!!!do i=1, n_snapshot
    !!!!!!    do j=1, Nmodes
    !!!!!!        tmp(j,i)=Fcoeff2(i,j)
    !!!!!!    enddo
    !!!!!!enddo
    !!!!!!allocate(re_u(m_dim,n_snapshot), re_v(m_dim,n_snapshot))
    !!!!!!call MM(m_dim, Nmodes, Nmodes, n_snapshot, u, tmp, re_u)
    !!!!!!call MM(m_dim, Nmodes, Nmodes, n_snapshot, v, tmp, re_v)
    !!!!!!deallocate(tmp)
    !!!!!!do i=1, n_snapshot
    !!!!!!    re_u(:,i)=re_u(:,i)+u_avg &
    !!!!!!        +gam(n_snapshot*1+i)*ConFM(:,1)+gam2(n_snapshot*1+i)*wall_v1(:,1)
    !!!!!!    re_v(:,i)=re_v(:,i)+v_avg &
    !!!!!!        +gam(n_snapshot*1+i)*ConFM(:,2)+gam2(n_snapshot*1+i)*wall_v1(:,2)
    !!!!!!enddo
    !!!!!!i1=0
    !!!!!!do k=1, 0+1
    !!!!!!    do i=1, n_snapshot
    !!!!!!        i1=i1+1
    !!!!!!        t_read=t_start2+dt*dble(i)
    !!!!!!        write(time,*) floor(t_read)
    !!!!!!        time=adjustl(time)
    !!!!!!    
    !!!!!!        open(1,file='re_construction1/re_'//time(1:len_trim(time))//'.dat')
    !!!!!!        write(1,49)
    !!!!!!        write(1,47) floor(t_read), m_dim, nel
    !!!!!!        do j=1, m_dim
    !!!!!!            write(1,48) x(j), y(j)+0.1d0*dcos(2.d0*pi*nd_freq1*dt_real1*dble(i)), &
    !!!!!!            !write(1,48) x(j), y(j)+0.1d0*dcos(2.d0*pi*0.1957*1.2d0*dt_real*dble(i)), &     
    !!!!!!                re_u(j,i1)*Uref, re_v(j,i1)*Uref, 0.d0
    !!!!!!        enddo
    !!!!!!        
    !!!!!!        do j=1, nel
    !!!!!!            write(1,*) conv(j,1), conv(j,2), conv(j,3)
    !!!!!!        enddo
    !!!!!!        close(1)
    !!!!!!    enddo
    !!!!!!enddo
    deallocate(re_u,re_v)
    
    Xpos=x(:); Ypos=y(:)
    AVGu=u_avg; AVGv=v_avg; AVGp=0.d0
    PODu=u(:,1:Nmodes); PODv=v(:,1:Nmodes); PODp=0.d0
    yval(:)=Fcoeff(1,:)
    allocate(CON_u(m_dim,1), CON_v(m_dim,1))
    CON_u(:,1)=ConFM(:,1); CON_v(:,1)=ConFM(:,2)
    deallocate( u, v, u1, v1, x, y, u_avg, v_avg, ConFM )
    
49	FORMAT('variables = "X","Y","u","v","p"')
47  FORMAT('zone t="',I9,'", N=  ',I9,' , E=  ',I9,' , DATAPACKING = POINT, ZONETYPE = FETRIANGLE')
!47	FORMAT('zone t="',I9,'", i=321, j=140 , DATAPACKING = POINT')
48  FORMAT(F23.16,F23.16,F23.16,F23.16,F23.16,F23.16 )	   
449	FORMAT(100F23.16) 
    end subroutine
    