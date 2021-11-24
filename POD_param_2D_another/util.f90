
    
    
    subroutine eig_mat(m,n,a,vt,s)
    implicit none
    
    integer m, n, i
    double precision a(m,n), vt(n,n), s(n), temp(n), temp_v(n,n)
    INTEGER          INFO, LWORK
    double precision, allocatable :: u(:,:), work(:)
    
    allocate( WORK( M*M ) )
    
    
    call DSYEV( 'V', 'L', N, A, N, S, WORK, M*M, INFO )
    
    ! n_snap_t = m,n_snap_t = n, tmp = a, u_vec= vt, u_f = s
    ! JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO
    ! N : order of matrix
    ! A : 
    
    !print*, info;pause;
    
    vt=A ! eigenvector
    
    
    
    !if (abs(s(n))>abs(s(1))) then ! Eigenvalue in the negative
    !    do i = 1,n
    !        temp(i) = s(n-i+1);
    !        temp_v(:,i) = VT(:,n-i+1);
    !    enddo
    !    do i = 1,n
    !        s(i) = temp(i);
    !        VT = temp_v   ! eigenvalue switched
    !    enddo
    !endif
    !print*, temp;pause;
    !print*, s;pause;
    
    !s=s;
    
    
    deallocate(work)
    end subroutine
    
    
    SUBROUTINE MM(N1, N2, M1, M2, MAT1, MAT2, NEW_MAT)
    IMPLICIT NONE
    INTEGER I, J, K
    
    INTEGER N1, N2, M1, M2
    DOUBLE PRECISION MAT1(N1, N2), MAT2(M1, M2)
    DOUBLE PRECISION NEW_MAT(N1,M2)
    NEW_MAT=0.D0
    IF (N2 .NE. M1) THEN
        WRITE(*,*) 'PRODUCT IS NOT POSSIBLE !!'
        PAUSE
    ENDIF

    DO I=1, N1
        DO K=1,M2
            DO J=1, N2
                NEW_MAT(I,K)=NEW_MAT(I,K)+MAT1(I,J)*MAT2(J,K)
                
                !print*, NEW_MAT(I,K)
                !print*, MAT1(I,J)
                !print*, MAT2(J,K)
                !pause;
            ENDDO
        ENDDO    
    ENDDO
    
    end subroutine
    
    
    