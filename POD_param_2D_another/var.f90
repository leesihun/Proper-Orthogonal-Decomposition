
    module var_POD
        !integer, parameter :: n_snapshot=400, n_snapshot2=0, Nmodes=8, Ncntrl=1
        !integer, parameter :: nnode=63388, nelem=126178
        !double precision, parameter :: Re=200.0
        
        integer, parameter :: n_snapshot=290   ! total number of snapshot = 2*n_snapshot
        integer, parameter :: Nmodes=13        ! number of POD modes
        integer, parameter :: nnode = 55634, nelem = 110581 ! Cylinder nnode=55634, nelem=110581
        double precision, parameter :: Re=1000.0
        
        integer, allocatable :: conv(:,:), ia(:,:)
        double precision, allocatable :: weight(:,:)
        double precision, allocatable :: PODu(:,:), PODv(:,:), PODp(:,:)
        double precision, allocatable :: AVGu(:), AVGv(:), AVGp(:)
        double precision, allocatable :: Xpos(:), Ypos(:)
        double precision, allocatable :: A(:), B(:,:), C(:,:,:), Ag(:), Bg(:,:)
        double precision, allocatable :: D(:,:), E(:), F(:,:), G(:,:,:)
        
        double precision, allocatable :: CON_u(:,:), CON_v(:,:), KMg(:)
        double precision, allocatable :: HM(:,:), JM(:,:,:), KM(:,:), LM(:,:,:),LMg(:,:)
        
        double precision, allocatable :: Gvel_u(:), Gvel_v(:)
        double precision yval(Nmodes), val_v
        integer, parameter :: n_snapshot2=0, Ncntrl=1   ! not used in this program
    end module
    
    
    
    
