Program Quasimagnetostaticsolver3D

    !LICENSE FILE is included in this folder or in the upper folder
    !Assume all SI units throughout
    !this program solves the 3D quasi-magnetostatic equations in the Coulomb gauge form
    !It uses both of these equations:
    !1: del cross (mu inverse times del cross A) + sigma*del(V) - gauge penalt term + sigma * del A/ del t - sigma * v cross del cross A = J_applied (zero here)
    !2: del dot J_total = 0 = del dot (sigma*del(V) + sigma*del A/ del t - sigma*v cross del cross A - J_applied (zero here))
    
    implicit none

    ! Precision definition for double precision (64-bit real)
    integer, parameter :: dp = selected_real_kind(15, 307)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Setup all constants and arrays !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !constants
    real(dp), parameter :: c = 299792458.0_dp                  !speed of light
    real(dp), parameter :: mu_0 = 1.25663706E-6_dp              !permeability of free space

    !basic parameters
    real(dp) :: del_t                                           !setup delta t
    integer :: time_steps                                       !number of time steps
    real(dp) :: del_x,del_y,del_z                               !step sizes in each direction x,y,z
    integer :: x_size,y_size,z_size                             !size of grid (# of cubes)
    integer :: Nx_nodes,Ny_nodes,Nz_nodes                       !size of grid (# of nodes)
    integer :: num_materials                                    !number of materials user submits

    !helper variables
    integer :: i,j,k,kk,rr, rr_rel, counter                     !used for looping
    real(dp) :: dAx_dy,dAx_dz                                   !temps for setting up H,B calculations
    real(dp) :: dAy_dx,dAy_dz                                   !temps for setting up H,B calculations
    real(dp) :: dAz_dx,dAz_dy                                   !temps for setting up H,B calculations
    integer :: materials_id                                     !temporarily tracks the material ID
    real(dp) :: x_pos,y_pos,z_pos                               !temp vars in solvers
    real(dp) :: w                                               !temp var in solvers
    integer :: tt                                               !used to mod time steps for saving fields every tt time step
    integer :: exact_nn, num_interior, num_boundary             !temp vars used in solvers
    integer :: nn_counter                                       !temp var used in solvers
    character(3) :: output_slice                                !output slice type x,y,z for output field slice
    character(3) :: save_all_fields                             !save all fields - 'yes' or 'no'
    integer :: location_slice                                   !slice location for output

    !time keeping variables
    integer :: clock_time_start                                 !for timing of simulation
    integer :: clock_time_end                                   !for timing of simulation
    integer :: clock_rate                                       !for timing of simulation

    !file variables
    character(len=50) :: input_file_name                        !name of file to read into this programn for sim parameters
    character(len=50) :: input_geom_name                        !name of file to read into this program for geometry

    !initial and time dependent boundary conditions
    real(dp) :: Bx_in,By_in,Bz_in                               !Background displacement field
    integer :: type_traj_x,type_traj_y,type_traj_z              !built in trajectory type
    real(dp) :: init_vx,init_vy,init_vz                         !initial velocity in x,y,z directions for traj type 1
    real(dp) :: const_ax,const_ay,const_az                      !acceleration in x,y,z direction for traj type 1
    integer :: cax_time,cay_time,caz_time                       !time steps for acceleration in x,y,z direction for traj type 1
    real(dp) :: solve_vx, solve_vy, solve_vz                    !velocities at a given time step
    real(dp) :: o_solve_vx, o_solve_vy, o_solve_vz              !velocities at a given previous time step (n-1) - for speed boost
    real(dp) :: sigx0,sigx1,sigy0,sigy1,sigz0,sigz1             !helper vars to clean up V section
    real(dp) :: temp_acc_x,temp_acc_y,temp_acc_z                !helper vars for stability purposes
    real(dp) :: gxx,gxy,gxz,gyx,gyy,gyz,gzx,gzy,gzz             !gradients for B fields, first index is field component, and then direction

    !Field arrays
    real(dp), allocatable :: Ax(:,:,:)                          !Magnetic vector potential in x direction
    real(dp), allocatable :: Ay(:,:,:)                          !Magnetic vector potential in y direction
    real(dp), allocatable :: Az(:,:,:)                          !Magnetic vector potential in z direction
    real(dp), allocatable :: V(:,:,:)                           !Electric scalar potential
    real(dp), allocatable :: Hx(:,:,:)                          !Magnetic field x component
    real(dp), allocatable :: Hy(:,:,:)                          !Magnetic field y component
    real(dp), allocatable :: Hz(:,:,:)                          !Magnetic field y component
    real(dp), allocatable :: Hx_out(:,:)                        !Magnetic field x component - only output part of the fields as a 2D slice
    real(dp), allocatable :: Hy_out(:,:)                        !Magnetic field y component - only output part of the fields as a 2D slice
    real(dp), allocatable :: Hz_out(:,:)                        !Magnetic field z component - only output part of the fields as a 2D slice

    !Material arrays - uses yee cell positioning like in FDTD
    real(dp), allocatable :: Mu_cells_x(:,:,:)                  !Permeability of the cells in x direction
    real(dp), allocatable :: Mu_cells_y(:,:,:)                  !Permeability of the cells in y direction
    real(dp), allocatable :: Mu_cells_z(:,:,:)                  !Permeability of the cells in z direction
    real(dp), allocatable :: Mu_x(:,:,:)                        !Permeability of the Yee edges in x direction
    real(dp), allocatable :: Mu_y(:,:,:)                        !Permeability of the Yee edges in y direction
    real(dp), allocatable :: Mu_z(:,:,:)                        !Permeability of the Yee edges in z direction
    real(dp), allocatable :: Sigma_cells_x(:,:,:)               !Electrical conductivity of the cells in x direction
    real(dp), allocatable :: Sigma_cells_y(:,:,:)               !Electrical conductivity of the cells in y direction
    real(dp), allocatable :: Sigma_cells_z(:,:,:)               !Electrical conductivity of the cells in z direction
    real(dp), allocatable :: Sigma_x(:,:,:)                     !Electrical conductivity of the Yee edgesa in x direction
    real(dp), allocatable :: Sigma_y(:,:,:)                     !Electrical conductivity of the Yee edgesa in y direction
    real(dp), allocatable :: Sigma_z(:,:,:)                     !Electrical conductivity of the Yee edgesa in z direction

    !Other arrays
    real(dp), allocatable :: Fixed_x(:,:,:)                     !Fixed edges by boundary conditions in the solver array (a_mat*x_vec=b_vec)
    real(dp), allocatable :: Fixed_y(:,:,:)                     !Fixed edges by boundary conditions in the solver array (a_mat*x_vec=b_vec)
    real(dp), allocatable :: Fixed_z(:,:,:)                     !Fixed edges by boundary conditions in the solver array (a_mat*x_vec=b_vec)
    real(dp), allocatable :: Fixed_v(:,:,:)                     !Fixed nodes by boundary conditions in the solver array (a_mat*x_vec=b_vec)

    real(dp), allocatable :: material_properties(:,:,:)         !Material properties data
    real(dp), allocatable :: input_geom_file(:,:,:)             !Array of integers (as floats) indicating geometry
    integer, allocatable :: temp_mat_id(:)                      !Used for geometry read in for identifying properties associated with a mat ID

    !PARDISO specific arrays
    !originally wanted the option to use different lengths for ja and a_mat for advective vs normal - but currently use identical lengths
    integer(8) :: pt(64)
    integer :: maxfct, mnum, mtype, phase, n_dims, n_dims_1vec, nrhs, error, msglvl
    integer :: iparm(64)
    integer, allocatable :: ia(:), ja(:)
    real(8), allocatable :: a_mat(:), b_vec(:), x_vec(:)       ! Matches real(dp) precision

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Inputs read in section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call system_clock(clock_time_start,clock_rate)
    write(*, '("Reading in the inputs file and setting up the simulation...")')

    call get_command_argument(1, input_file_name)
    open(1,file=trim(input_file_name),status="old",action="read")

        read(1,*) x_size,y_size,z_size
        read(1,*) del_x,del_y,del_z
        read(1,*) del_t
        read(1,*) time_steps
        read(1,*) tt
        read(1,*) output_slice,location_slice
        read(1,*) save_all_fields
        read(1,*) Bx_in,By_in,Bz_in
        read(1,*) gxx,gxy,gxz
        read(1,*) gyx,gyy,gyz
        read(1,*) gzx,gzy,gzz
        read(1,*) type_traj_x
        if (type_traj_x==1) then
            read(1,*) init_vx,const_ax,cax_time
        else
            STOP "ERROR: TYPE OF TRAJECTORY IN X DIRECTION IS NOT SUPPORTED"
        end if
        read(1,*) type_traj_y
        if (type_traj_y==1) then
            read(1,*) init_vy,const_ay,cay_time
        else
            STOP "ERROR: TYPE OF TRAJECTORY IN Y DIRECTION IS NOT SUPPORTED"
        end if
        read(1,*) type_traj_z
        if (type_traj_z==1) then
            read(1,*) init_vz,const_az,caz_time
        else
            STOP "ERROR: TYPE OF TRAJECTORY IN Z DIRECTION IS NOT SUPPORTED"
        end if

        read(1,*) num_materials
        if (num_materials>0) then
            allocate(material_properties(num_materials,3,3))
        end if
        do i=1, num_materials
            read(1,*) material_properties(i,1,1)
            read(1,*) material_properties(i,2,1), material_properties(i,2,2), material_properties(i,2,3)
            read(1,*) material_properties(i,3,1), material_properties(i,3,2), material_properties(i,3,3)
        end do

        read(1,*) input_geom_name

    close(1)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Assign some variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Nx_nodes = x_size + 1
    Ny_nodes = y_size + 1
    Nz_nodes = z_size + 1
    n_dims = Nx_nodes * Ny_nodes * Nz_nodes * 4 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes - Nx_nodes*Ny_nodes
    n_dims_1vec = Nx_nodes * Ny_nodes * Nz_nodes

    !exact assumes only edges are fixed - it can be too large at cost of memory, just not too small wll throw an error statement
    num_boundary = (2 * Nx_nodes + 2 * Ny_nodes - 4) * Nz_nodes   !2x walls + 2y walls -4 for double counting corners times the number of Z layers there are
    num_interior = ((Nx_nodes - 2) * (Ny_nodes - 2)) * Nz_nodes   !remove 2 from the distances and find the area times number of Z layers there are
    exact_nn = (num_boundary * 1) + (num_interior * (15*3 + 37))  !add the numbers together - boundary has 1 entry and num_interior has 15*3 + 37 entries

    !intialize o_solves - not needed in steady state solver - used to boost speed in main time loop if speed didn't change
    o_solve_vx = trajectory_1(init_vx, const_ax, cax_time, 0, del_t)
    o_solve_vy = trajectory_1(init_vy, const_ay, cay_time, 0, del_t)
    o_solve_vz = trajectory_1(init_vz, const_az, caz_time, 0, del_t)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Setting up the geometry !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !allocate material arrays and intialize them as is appropriate - both cells and Yee grid positionings

    allocate(Mu_cells_x(x_size,y_size,z_size))
    allocate(Mu_cells_y(x_size,y_size,z_size))
    allocate(Mu_cells_z(x_size,y_size,z_size))

    allocate(Mu_x(Nx_nodes,Ny_nodes,Nz_nodes))
    allocate(Mu_y(Nx_nodes,Ny_nodes,Nz_nodes))
    allocate(Mu_z(Nx_nodes,Ny_nodes,Nz_nodes))

    allocate(Sigma_cells_x(x_size,y_size,z_size))
    allocate(Sigma_cells_y(x_size,y_size,z_size))
    allocate(Sigma_cells_z(x_size,y_size,z_size))

    allocate(Sigma_x(Nx_nodes,Ny_nodes,Nz_nodes))
    allocate(Sigma_y(Nx_nodes,Ny_nodes,Nz_nodes))
    allocate(Sigma_z(Nx_nodes,Ny_nodes,Nz_nodes))

    Mu_cells_x(:,:,:)=mu_0
    Mu_cells_y(:,:,:)=mu_0
    Mu_cells_z(:,:,:)=mu_0

    Mu_x(:,:,:)=mu_0
    Mu_y(:,:,:)=mu_0
    Mu_z(:,:,:)=mu_0

    Sigma_cells_x(:,:,:)=0.0
    Sigma_cells_y(:,:,:)=0.0
    Sigma_cells_z(:,:,:)=0.0

    Sigma_x(:,:,:)=0.0
    Sigma_y(:,:,:)=0.0
    Sigma_z(:,:,:)=0.0

    !now read in the geometry file and assign the cells the corresponding materials

    allocate(input_geom_file(x_size,y_size,z_size))       
    open(unit=2, file=trim(input_geom_name), access='stream', status='old')
        read(2) input_geom_file
    close(2)
    
    allocate(temp_mat_id(0:maxval(int(material_properties(:, 1, 1)))))
    temp_mat_id(:) = 0
    do rr=1, num_materials
        temp_mat_id(int(material_properties(rr, 1, 1))) = rr
    end do
    
    do k=1, z_size
        do j=1, y_size
            do i=1, x_size
                materials_id = temp_mat_id(int(input_geom_file(i,j,k)))
                if (materials_id /= 0) then
                    Mu_cells_x(i,j,k)=material_properties(materials_id,2,1)*mu_0
                    Mu_cells_y(i,j,k)=material_properties(materials_id,2,2)*mu_0
                    Mu_cells_z(i,j,k)=material_properties(materials_id,2,3)*mu_0
                    Sigma_cells_x(i,j,k)=material_properties(materials_id,3,1)
                    Sigma_cells_y(i,j,k)=material_properties(materials_id,3,2)
                    Sigma_cells_z(i,j,k)=material_properties(materials_id,3,3)
                end if
            end do
        end do
    end do

    !now determine the Yee positioning for sigma and mu, respectively.
    !Don't worry about edges because they need to be vacuum properties anyway as was defaulted above
    do k=2, z_size
        do j=2, y_size
            do i=2, x_size
                Mu_x(i,j,k)=2.0*Mu_cells_x(i,j,k)*Mu_cells_x(i-1,j,k)/(Mu_cells_x(i,j,k)+Mu_cells_x(i-1,j,k))
                Mu_y(i,j,k)=2.0*Mu_cells_y(i,j,k)*Mu_cells_y(i,j-1,k)/(Mu_cells_y(i,j,k)+Mu_cells_y(i,j-1,k))
                Mu_z(i,j,k)=2.0*Mu_cells_z(i,j,k)*Mu_cells_z(i,j,k-1)/(Mu_cells_z(i,j,k)+Mu_cells_z(i,j,k-1))
                Sigma_x(i,j,k)=(Sigma_cells_x(i,j,k)+Sigma_cells_x(i,j,k-1)+Sigma_cells_x(i,j-1,k)+Sigma_cells_x(i,j-1,k-1))/4.0
                Sigma_y(i,j,k)=(Sigma_cells_y(i,j,k)+Sigma_cells_y(i-1,j,k)+Sigma_cells_y(i,j,k-1)+Sigma_cells_y(i-1,j,k-1))/4.0
                Sigma_z(i,j,k)=(Sigma_cells_z(i,j,k)+Sigma_cells_z(i-1,j,k)+Sigma_cells_z(i,j-1,k)+Sigma_cells_z(i-1,j-1,k))/4.0
            end do
        end do
    end do

    !deallocate some memory we don't need anymore
    deallocate(temp_mat_id)
    deallocate(input_geom_file)    
    deallocate(Sigma_cells_x)
    deallocate(Sigma_cells_y)
    deallocate(Sigma_cells_z)
    deallocate(Mu_cells_x)
    deallocate(Mu_cells_y)
    deallocate(Mu_cells_z)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Setup fields and solver settings !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !allocation and intialization

    !vector potential shares the same positions as E fields do in the Yee grid - allocated node size though for convenience
    allocate(Ax(Nx_nodes,Ny_nodes,Nz_nodes))
    allocate(Ay(Nx_nodes,Ny_nodes,Nz_nodes))
    allocate(Az(Nx_nodes,Ny_nodes,Nz_nodes))

    !scalar potential - sits at nodes
    allocate(V(Nx_nodes,Ny_nodes,Nz_nodes))
    
    !normal Yee grid positions for H fields - allocated node size though for convenience
    allocate(Hx(Nx_nodes,Ny_nodes,Nz_nodes))
    allocate(Hy(Nx_nodes,Ny_nodes,Nz_nodes))
    allocate(Hz(Nx_nodes,Ny_nodes,Nz_nodes))
    if (TRIM(output_slice)=='x') then
        allocate(Hx_out(y_size,z_size))
        allocate(Hy_out(y_size,z_size))
        allocate(Hz_out(y_size,z_size))
    end if
    if (TRIM(output_slice)=='y') then
        allocate(Hx_out(x_size,z_size))
        allocate(Hy_out(x_size,z_size))
        allocate(Hz_out(x_size,z_size))
    end if
    if (TRIM(output_slice)=='z') then
        allocate(Hx_out(x_size,y_size))
        allocate(Hy_out(x_size,y_size))
        allocate(Hz_out(x_size,y_size))
    end if

    !Each node as the potential to be fixed, though we will only allow to fix the edges (or nodes for V) in this solver for now
    allocate(Fixed_x(x_size,Ny_nodes,Nz_nodes))
    allocate(Fixed_y(Nx_nodes,y_size,Nz_nodes))
    allocate(Fixed_z(Nx_nodes,Ny_nodes,z_size))
    allocate(Fixed_v(Nx_nodes,Ny_nodes,Nz_nodes))

    !solver specifc array allocations
    allocate(ia(n_dims+1))
    allocate(b_vec(n_dims))
    allocate(x_vec(n_dims))
    allocate(ja(exact_nn))
    allocate(a_mat(exact_nn))

    !array initializations

    !vector fields
    Ax(:,:,:)=0.0
    Ay(:,:,:)=0.0
    Az(:,:,:)=0.0

    !scalar fields
    V(:,:,:)=0.0

    !H fields
    Hx(:,:,:)=0.0
    Hy(:,:,:)=0.0
    Hz(:,:,:)=0.0
    Hx_out(:,:)=0.0
    Hy_out(:,:)=0.0
    Hz_out(:,:)=0.0

    !bcs
    Fixed_x(:,:,:)=0.0
    Fixed_y(:,:,:)=0.0
    Fixed_z(:,:,:)=0.0
    Fixed_v(:,:,:)=0.0

    !solver specific
    b_vec(:)=0.0
    ia(:)=0.0
    ja(:)=0.0
    a_mat(:)=0.0

    !intialize the outside edges to have fixed values - will reference this array to know if it was fixed or not
    !Only outer edge is required to be fixed but we could fix as many as we want moving inwards - I'm doing entire cells right now
    do k = 1, Nz_nodes
        do j = 1, Ny_nodes
            do i = 1, x_size
                if (j <= 2 .or. j >= Ny_nodes-1 .or. k <= 2 .or. k >= Nz_nodes-1 .or. i==1 .or. i==x_size) then
                    Fixed_x(i,j,k) = 1.0
                end if
            end do
        end do
    end do
    do k = 1, Nz_nodes
        do j = 1, y_size
            do i = 1, Nx_nodes
                if (i <= 2 .or. i >= Nx_nodes-1 .or. k <= 2 .or. k >= Nz_nodes-1 .or. j==1 .or. j==y_size) then
                    Fixed_y(i,j,k) = 1.0
                end if
            end do
        end do
    end do
    do k = 1, z_size
        do j = 1, Ny_nodes
            do i = 1, Nx_nodes
                if (i <= 2 .or. i >= Nx_nodes-1 .or. j <= 2 .or. j >= Ny_nodes-1 .or. k==1 .or. k==z_size) then
                    Fixed_z(i,j,k) = 1.0
                end if
            end do
        end do
    end do
    do k = 1, Nz_nodes
        do j = 1, Ny_nodes
            do i = 1, Nx_nodes
                if (i <= 2 .or. i >= Nx_nodes-1 .or. j <= 2 .or. j >= Ny_nodes-1 .or. k <= 2 .or. k >= Nz_nodes-1) then
                    Fixed_v(i,j,k) = 1.0
                end if
            end do
        end do
    end do

    !PARDISO Configuration Settings
    pt = 0                ! Initialize 64-bit handle markers directly as zeros
    iparm = 0
    iparm(1) = 1          ! Use customized configurations
    iparm(2) = 2          ! Nested dissection fill-in reducing ordering
    iparm(28) = 0         ! Double precision execution
    mtype = 11            ! Real structurally unsymmetric matrix
    maxfct = 1; mnum = 1; nrhs = 1; msglvl = 0

    !To save RAM, we will save at each time step division chosen by user to the binary - this will avoid saving Hx,Hy,Hz in time
    open(3, file="Hx.bin", form="unformatted",action="write",status="replace")
    open(4, file="Hy.bin", form="unformatted",action="write",status="replace")
    open(5, file="Hz.bin", form="unformatted",action="write",status="replace")

    !if saving all fields
    if (TRIM(save_all_fields)=='yes') then
        open(6, file="Hx_all.bin", form="unformatted",action="write",status="replace")
        open(7, file="Hy_all.bin", form="unformatted",action="write",status="replace")
        open(8, file="Hz_all.bin", form="unformatted",action="write",status="replace")
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Initialize via steady state (v=0 static or v/=0 advection-diffusion solver) !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call system_clock(clock_time_end, clock_rate)
    write(*, '("Setup completed in ", F6.2, " seconds.")') &
        real(clock_time_end - clock_time_start) / real(clock_rate)
    write(*, '("Intializing steady state...")')
    call system_clock(clock_time_start)    

    !determine trajectory (velocity) at time=0
    solve_vx = trajectory_1(init_vx, const_ax, cax_time, 0, del_t)
    solve_vy = trajectory_1(init_vy, const_ay, cay_time, 0, del_t)
    solve_vz = trajectory_1(init_vz, const_az, caz_time, 0, del_t)

    !Setting up b and CSR structures (a_mat_adv, ja, ia) for initialization procedure
    nn_counter = 1

    !Loop sequentially through rows from 1 to n_dims so ia is built perfectly
    !the rows are all of the equations - there is an equation for each Ax,Ay,Az at each position
    !Easiest to loop through all Ax, then all Ay, then all Az

    !steady state has del A/ del t set to zero for steady state purposes
    !full paridso phase 13 is used for setup and solve

    !!!!!!!!!!!!!!!!!!!!
    !!!!First is Ax!!!!!
    !!!!!!!!!!!!!!!!!!!!
    do rr = 1, n_dims_1vec-Ny_nodes*Nz_nodes
        
        !tell the solver where each row starts
        !nn_counter gets udpated each time a value is added to a_mat below, always starts at 1 as set above
        ia(rr) = nn_counter

        !Ordering here is very important for PARDISO, we need to add them in order down the row
        i = mod(rr - 1, x_size) + 1
        j = mod((rr - 1) / x_size, Ny_nodes) + 1
        k = (rr - 1) / (Ny_nodes * x_size) + 1

        !First check if it's a fixed entry (boundary condition) and if so we know and set A at that position, then move on (cycle)
        if (Fixed_x(i, j, k) == 1.0) then
            x_pos = (i - 1) * del_x
            y_pos = (j - 1) * del_y
            z_pos = (k - 1) * del_z
            b_vec(rr) = 0.5 * By_in * (z_pos) - 0.5 * Bz_in * (y_pos) + &
            1/3.0 * ((gyx*x_pos + gyy*y_pos + gyz*z_pos)*z_pos - (gzx*x_pos + gzy*y_pos + gzz*z_pos)*y_pos)
            a_mat(nn_counter) = 1.0
            ja(nn_counter) = rr
            nn_counter = nn_counter + 1
            cycle
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!First Ax self entries!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Ax(i,j,k-1)
        w = -1.0 / Mu_y(i,j,k-1) / (del_z ** 2)  + (sigma_x(i,j,k) * solve_vz) / (2.0 * del_z)
        kk = i + (j - 1) * x_size + (k - 2) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i,j-1,k)
        w = -1.0 / Mu_z(i,j-1,k) / (del_y ** 2) + (sigma_x(i,j,k) * solve_vy) / (2.0 * del_y)
        kk = i + (j - 2) * x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !GAUGE PENALTY: Ax(i-1,j,k) - Added to enforce -d/dx( 1/mu * dAx/dx )
        w = -1.0 / Mu_x(i-1,j,k) / (del_x ** 2)
        kk = (i - 1) + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i,j,k) (Diagonal spatial + gauge penalty)
        w = (1.0 / Mu_z(i,j,k) + 1.0 / Mu_z(i,j-1,k)) / (del_y ** 2) &
          + (1.0 / Mu_y(i,j,k) + 1.0 / Mu_y(i,j,k-1)) / (del_z ** 2) &
          + (1.0 / Mu_x(i,j,k) + 1.0 / Mu_x(i-1,j,k)) / (del_x ** 2)
        kk = i + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !GAUGE PENALTY: Ax(i+1,j,k) - Added to enforce -d/dx( 1/mu * dAx/dx )
        w = -1.0 / Mu_x(i,j,k) / (del_x ** 2)
        kk = (i + 1) + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i,j+1,k)
        w = -1.0 / Mu_z(i,j,k) / (del_y ** 2) - (sigma_x(i,j,k) * solve_vy) / (2.0 * del_y)
        kk = i + (j + 0) * x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i,j,k+1)
        w = -1.0 / Mu_y(i,j,k) / (del_z ** 2) - (sigma_x(i,j,k) * solve_vz) / (2.0 * del_z)
        kk = i + (j - 1) * x_size + (k + 0) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!Second Ax has Ay entries
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Ay(i,j-1,k)
        w = 1.0 / Mu_z(i,j-1,k) / (del_x *del_y) - (sigma_x(i,j,k) * solve_vy) / (2.0 * del_x)
        kk = i + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec-Ny_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i+1,j-1,k)
        w = -1.0 / Mu_z(i,j-1,k) / (del_x *del_y) + (sigma_x(i,j,k) * solve_vy) / (2.0 * del_x)
        kk = i + 1 + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec-Ny_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i,j,k)
        w = -1.0 / Mu_z(i,j,k) / (del_x *del_y) - (sigma_x(i,j,k) * solve_vy) / (2.0 * del_x)
        kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec-Ny_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i+1,j,k)
        w = 1.0 / Mu_z(i,j,k) / (del_x *del_y) + (sigma_x(i,j,k) * solve_vy) / (2.0 * del_x)
        kk = i + 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec-Ny_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!Third Ax has Az entries
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Az(i,j,k-1)
        w = 1.0 / Mu_y(i,j,k-1) / (del_x *del_z) - (sigma_x(i,j,k) * solve_vz) / (2.0 * del_x)
        kk = i + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + n_dims_1vec*2-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i+1,j,k-1)
        w = -1.0 / Mu_y(i,j,k-1) / (del_x *del_z) + (sigma_x(i,j,k) * solve_vz) / (2.0 * del_x)
        kk = i + 1 + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + n_dims_1vec*2-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i,j,k)
        w = -1.0 / Mu_y(i,j,k) / (del_x *del_z) - (sigma_x(i,j,k) * solve_vz) / (2.0 * del_x)
        kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i+1,j,k)
        w = 1.0 / Mu_y(i,j,k) / (del_x *del_z) + (sigma_x(i,j,k) * solve_vz) / (2.0 * del_x)
        kk = i + 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!Lastly Ax has V entries
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !V(i,j,k)
        w = -Sigma_x(i,j,k) / del_x
        kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !V(i+1,j,k)
        w = Sigma_x(i,j,k) / del_x
        kk = i + 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Safety check in case I miscalculated the most possible entries
        if (nn_counter > exact_nn) then
            write(*,*) "Error: a_mat allocation exceeded! Increase multiplier."
            stop
        end if

    end do

    !!!!!!!!!!!!!!!!!!!!
    !!!!Second is Ay!!!!
    !!!!!!!!!!!!!!!!!!!!
    do rr = n_dims_1vec - Ny_nodes*Nz_nodes + 1, n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes

        !tell the solver where each row starts
        !nn_counter gets udpated each time a value is added to a_mat below, always starts at 1 as set above
        ia(rr) = nn_counter

        !Ordering here is very important for PARDISO, we need to add them in order down the row

        rr_rel = rr - (n_dims_1vec - Ny_nodes*Nz_nodes)
        i = mod(rr_rel - 1, Nx_nodes) + 1
        j = mod((rr_rel - 1) / Nx_nodes, y_size) + 1
        k = (rr_rel - 1) / (y_size * Nx_nodes) + 1

        !First check if it's a fixed entry (boundary condition) and if so we know and set A at that position, then move on (cycle)
        if (Fixed_y(i, j, k) == 1.0) then
            x_pos = (i - 1) * del_x
            y_pos = (j - 1) * del_y
            z_pos = (k - 1) * del_z
            b_vec(rr) = 0.5 * Bz_in * (x_pos) - 0.5 * Bx_in * (z_pos) + &
            1/3.0 * ((gzx*x_pos + gzy*y_pos + gzz*z_pos)*x_pos - (gxx*x_pos + gxy*y_pos + gxz*z_pos)*z_pos)
            a_mat(nn_counter) = 1.0
            ja(nn_counter) = rr
            nn_counter = nn_counter + 1
            cycle
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!First Ay has Ax entries!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Ax(i-1,j,k)
        w = 1.0 / Mu_z(i-1,j,k) / (del_x * del_y) - (Sigma_y(i,j,k) * solve_vx) / (2.0 * del_y)
        kk = i - 1 + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i,j,k)
        w = -1.0 / Mu_z(i,j,k) / (del_x * del_y) - (Sigma_y(i,j,k) * solve_vx) / (2.0 * del_y)
        kk = i + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i-1,j+1,k)
        w = -1.0 / Mu_z(i-1,j,k) / (del_x * del_y) + (Sigma_y(i,j,k) * solve_vx) / (2.0 * del_y)
        kk = i - 1 + (j + 0) * x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i,j+1,k)
        w = 1.0 / Mu_z(i,j,k) / (del_x * del_y) + (Sigma_y(i,j,k) * solve_vx) / (2.0 * del_y)
        kk = i + (j + 0) * x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!Second Ay self entries!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Ay(i,j,k-1)
        w = -1.0 / Mu_x(i,j,k-1) / (del_z ** 2) + (Sigma_y(i,j,k) * solve_vz) / (2.0 * del_z)
        kk = i + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * y_size + n_dims_1vec - Ny_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !GAUGE PENALTY: Ay(i,j-1,k) - Added to enforce -d/dy( 1/mu * dAy/dy )
        w = -1.0 / Mu_y(i,j-1,k) / (del_y ** 2)
        kk = i + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec - Ny_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i-1,j,k)
        w = -1.0 / Mu_z(i-1,j,k) / (del_x ** 2) + (Sigma_y(i,j,k) * solve_vx) / (2.0 * del_x)
        kk = i - 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec - Ny_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i,j,k) (Diagonal spatial + gauge penalty)
        w = (1.0 / Mu_z(i,j,k) + 1.0 / Mu_z(i-1,j,k)) / (del_x ** 2) &
          + (1.0 / Mu_x(i,j,k) + 1.0 / Mu_x(i,j,k-1)) / (del_z ** 2) &
          + (1.0 / Mu_y(i,j,k) + 1.0 / Mu_y(i,j-1,k)) / (del_y ** 2)
        kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec - Ny_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i+1,j,k)
        w = -1.0 / Mu_z(i,j,k) / (del_x ** 2) - (Sigma_y(i,j,k) * solve_vx) / (2.0 * del_x)
        kk = i + 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec - Ny_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !GAUGE PENALTY: Ay(i,j+1,k) - Added to enforce -d/dy( 1/mu * dAy/dy )
        w = -1.0 / Mu_y(i,j,k) / (del_y ** 2)
        kk = i + (j + 0) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec - Ny_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i,j,k+1)
        w = -1.0 / Mu_x(i,j,k) / (del_z ** 2) - (Sigma_y(i,j,k) * solve_vz) / (2.0 * del_z)
        kk = i + (j - 1) * Nx_nodes + (k + 0) * Nx_nodes * y_size + n_dims_1vec - Ny_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!Third Ay has Az entries!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Az(i,j,k-1)
        w = 1.0 / Mu_x(i,j,k-1) / (del_y * del_z) - (Sigma_y(i,j,k) * solve_vz) / (2.0 * del_y)
        kk = i + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + n_dims_1vec*2-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i,j+1,k-1)
        w = -1.0 / Mu_x(i,j,k-1) / (del_y * del_z) + (Sigma_y(i,j,k) * solve_vz) / (2.0 * del_y)
        kk = i + (j + 0) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + n_dims_1vec*2-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i,j,k)
        w = -1.0 / Mu_x(i,j,k) / (del_y * del_z) - (Sigma_y(i,j,k) * solve_vz) / (2.0 * del_y)
        kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i,j+1,k)
        w = 1.0 / Mu_x(i,j,k) / (del_y * del_z) + (Sigma_y(i,j,k) * solve_vz) / (2.0 * del_y)
        kk = i + (j + 0) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!Lastly Ay has V entries
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !V(i,j,k)
        w = -Sigma_y(i,j,k) / del_y
        kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !V(i,j+1,k)
        w = Sigma_y(i,j,k) / del_y
        kk = i + (j + 0) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Safety check in case I miscalculated the most possible entries
        if (nn_counter > exact_nn) then
            write(*,*) "Error: a_mat allocation exceeded! Increase multiplier."
            stop
        end if

    end do

    !!!!!!!!!!!!!!!!!!!!
    !!!Third is Az!!!!!!
    !!!!!!!!!!!!!!!!!!!!
    do rr = n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes + 1, n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
        
        !tell the solver where each row starts
        !nn_counter gets udpated each time a value is added to a_mat below, always starts at 1 as set above
        ia(rr) = nn_counter

        !Ordering here is very important for PARDISO, we need to add them in order down the row

        rr_rel = rr - (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
        i = mod(rr_rel - 1, Nx_nodes) + 1
        j = mod((rr_rel - 1) / Nx_nodes, Ny_nodes) + 1
        k = (rr_rel - 1) / (Ny_nodes * Nx_nodes) + 1

        !First check if it's a fixed entry (boundary condition) and if so we know and set A at that position, then move on (cycle)
        if (Fixed_z(i, j, k) == 1.0) then
            x_pos = (i - 1) * del_x
            y_pos = (j - 1) * del_y
            z_pos = (k - 1) * del_z
            b_vec(rr) = 0.5 * Bx_in * y_pos - 0.5 * By_in * x_pos + &
            1/3.0 * ((gxx*x_pos + gxy*y_pos + gxz*z_pos)*y_pos - (gyx*x_pos + gyy*y_pos + gyz*z_pos)*x_pos)
            a_mat(nn_counter) = 1.0
            ja(nn_counter) = rr
            nn_counter = nn_counter + 1
            cycle
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!First Az has Ax entries!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Ax(i-1,j,k)
        w = 1.0 / Mu_y(i-1,j,k) / (del_x * del_z) - (Sigma_z(i,j,k) * solve_vx) / (2.0 * del_z)
        kk = i - 1 + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i,j,k)
        w = -1.0 / Mu_y(i,j,k) / (del_x * del_z) - (Sigma_z(i,j,k) * solve_vx) / (2.0 * del_z)
        kk = i + (j - 1) *x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i-1,j,k+1)
        w = -1.0 / Mu_y(i-1,j,k) / (del_x * del_z) + (Sigma_z(i,j,k) * solve_vx) / (2.0 * del_z)
        kk = i - 1 + (j - 1) * x_size + (k + 0) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i,j,k+1)
        w = 1.0 / Mu_y(i,j,k) / (del_x * del_z) + (Sigma_z(i,j,k) * solve_vx) / (2.0 * del_z)
        kk = i + (j - 1) * x_size + (k + 0) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!Second Az has Ay entries
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Ay(i,j-1,k)
        w = 1.0 / Mu_x(i,j-1,k) / (del_y * del_z) - (Sigma_z(i,j,k) * solve_vy) / (2.0 * del_z)
        kk = i + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec- Ny_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i,j,k)
        w = -1.0 / Mu_x(i,j,k) / (del_y * del_z) - (Sigma_z(i,j,k) * solve_vy) / (2.0 * del_z)
        kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size +n_dims_1vec- Ny_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i,j-1,k+1)
        w = -1.0 / Mu_x(i,j-1,k) / (del_y * del_z) + (Sigma_z(i,j,k) * solve_vy) / (2.0 * del_z)
        kk = i + (j - 2) * Nx_nodes + (k + 0) * Nx_nodes * y_size + n_dims_1vec- Ny_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i,j,k+1)
        w = 1.0 / Mu_x(i,j,k) / (del_y * del_z) + (Sigma_z(i,j,k) * solve_vy) / (2.0 * del_z)
        kk = i + (j - 1) * Nx_nodes + (k + 0) * Nx_nodes * y_size + n_dims_1vec- Ny_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!Third Az self entries!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !GAUGE PENALTY: Az(i,j,k-1) - Added to enforce -d/dz( 1/mu * dAz/dz )
        w = -1.0 / Mu_z(i,j,k-1) / (del_z ** 2)
        kk = i + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i,j-1,k)
        w = -1.0 / Mu_x(i,j-1,k) / (del_y ** 2) + (Sigma_z(i,j,k) * solve_vy) / (2.0 * del_y)
        kk = i + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i-1,j,k)
        w = -1.0 / Mu_y(i-1,j,k) / (del_x ** 2) + (Sigma_z(i,j,k) * solve_vx) / (2.0 * del_x)
        kk = i - 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i,j,k) (Diagonal spatial + gauge penalty)
        w = (1.0 / Mu_y(i,j,k) + 1.0 / Mu_y(i-1,j,k)) / (del_x ** 2) &
          + (1.0 / Mu_x(i,j,k) + 1.0 / Mu_x(i,j-1,k)) / (del_y ** 2) &
          + (1.0 / Mu_z(i,j,k) + 1.0 / Mu_z(i,j,k-1)) / (del_z ** 2)
        kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i+1,j,k)
        w = -1.0 / Mu_y(i,j,k) / (del_x ** 2) -(Sigma_z(i,j,k) * solve_vx) / (2.0 * del_x)
        kk = i + 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i,j+1,k)
        w = -1.0 / Mu_x(i,j,k) / (del_y ** 2) - (Sigma_z(i,j,k) * solve_vy) / (2.0 * del_y)
        kk = i + (j + 0) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !GAUGE PENALTY: Az(i,j,k+1) - Added to enforce -d/dz( 1/mu * dAz/dz )
        w = -1.0 / Mu_z(i,j,k) / (del_z ** 2)
        kk = i + (j - 1) * Nx_nodes + (k + 0) * Nx_nodes * Ny_nodes + n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!Lastly Az has V entries
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !V(i,j,k)
        w = -Sigma_z(i,j,k) / del_z
        kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !V(i,j,k+1)
        w = Sigma_z(i,j,k) / del_z
        kk = i + (j - 1) * Nx_nodes + (k + 0) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Safety check in case I miscalculated the most possible entries
        if (nn_counter > exact_nn) then
            write(*,*) "Error: a_mat allocation exceeded! Increase multiplier."
            stop
        end if

    end do

    !!!!!!!!!!!!!!!!!!!!
    !!!last is V!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!
    do rr = + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes + 1, n_dims
        
        !tell the solver where each row starts
        !nn_counter gets udpated each time a value is added to a_mat below, always starts at 1 as set above
        ia(rr) = nn_counter

        !Ordering here is very important for PARDISO, we need to add them in order down the row

        rr_rel = rr - (n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes)
        i = mod(rr_rel - 1, Nx_nodes) + 1
        j = mod((rr_rel - 1) / Nx_nodes, Ny_nodes) + 1
        k = (rr_rel - 1) / (Ny_nodes * Nx_nodes) + 1

        !Fix outer nodes to zero for convenience
        if (Fixed_V(i, j, k) == 1.0) then
            x_pos = (i - 1) * del_x
            y_pos = (j - 1) * del_y
            z_pos = (k - 1) * del_z
            b_vec(rr) = - 1.0 * (solve_vy * Bz_in - solve_vz * By_in) * x_pos &
            - 1.0 * (solve_vz * Bx_in - solve_vx * Bz_in) * y_pos &
            - 1.0 * (solve_vx * By_in - solve_vy * Bx_in) * z_pos &
            - 0.5 * (solve_vy * (gzx*x_pos + gzy*y_pos + gzz*z_pos) - solve_vz * (gyx*x_pos + gyy*y_pos + gyz*z_pos)) * x_pos &
            - 0.5 * (solve_vz * (gxx*x_pos + gxy*y_pos + gxz*z_pos) - solve_vx * (gzx*x_pos + gzy*y_pos + gzz*z_pos)) * y_pos &
            - 0.5 * (solve_vx * (gyx*x_pos + gyy*y_pos + gyz*z_pos) - solve_vy * (gxx*x_pos + gxy*y_pos + gxz*z_pos)) * z_pos
            a_mat(nn_counter) = 1.0
            ja(nn_counter) = rr
            nn_counter = nn_counter + 1
            cycle
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! div( sigma * (v x curl A) ) term added to the V equation !!!
        !!! (fully merged across x/y/z directions, sorted by kk)      !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        sigx0 = Sigma_x(i,j,k)
        sigx1 = Sigma_x(i-1,j,k)
        sigy0 = Sigma_y(i,j,k)
        sigy1 = Sigma_y(i,j-1,k)
        sigz0 = Sigma_z(i,j,k)
        sigz1 = Sigma_z(i,j,k-1)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! Ax cross entries !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Ax(i-1,j,k-1)
        w = sigx1 * solve_vz / (2.0 * del_z * del_x) - sigz1 * solve_vx / (2.0 * del_z ** 2)
        kk = (i - 1) + (j - 1) * x_size + (k - 2) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i,j,k-1)
        w = -sigx0 * solve_vz / (2.0 * del_z * del_x) - sigz1 * solve_vx / (2.0 * del_z ** 2)
        kk = i + (j - 1) * x_size + (k - 2) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i-1,j-1,k)
        w = sigx1 * solve_vy / (2.0 * del_y * del_x) - sigy1 * solve_vx / (2.0 * del_y ** 2)
        kk = (i - 1) + (j - 2) * x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i,j-1,k)
        w = -sigx0 * solve_vy / (2.0 * del_y * del_x) - sigy1 * solve_vx / (2.0 * del_y ** 2)
        kk = i + (j - 2) * x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i-1,j,k)
        w = (sigy0 + sigy1) * solve_vx / (2.0 * del_y ** 2) + (sigz0 + sigz1) * solve_vx / (2.0 * del_z ** 2)
        kk = (i - 1) + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i,j,k)
        w = (sigy0 + sigy1) * solve_vx / (2.0 * del_y ** 2) + (sigz0 + sigz1) * solve_vx / (2.0 * del_z ** 2)
        kk = i + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i-1,j+1,k)
        w = -sigx1 * solve_vy / (2.0 * del_y * del_x) - sigy0 * solve_vx / (2.0 * del_y ** 2)
        kk = (i - 1) + (j + 0) * x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i,j+1,k)
        w = sigx0 * solve_vy / (2.0 * del_y * del_x) - sigy0 * solve_vx / (2.0 * del_y ** 2)
        kk = i + (j + 0) * x_size + (k - 1) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i-1,j,k+1)
        w = -sigx1 * solve_vz / (2.0 * del_z * del_x) - sigz0 * solve_vx / (2.0 * del_z ** 2)
        kk = (i - 1) + (j - 1) * x_size + (k + 0) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ax(i,j,k+1)
        w = sigx0 * solve_vz / (2.0 * del_z * del_x) - sigz0 * solve_vx / (2.0 * del_z ** 2)
        kk = i + (j - 1) * x_size + (k + 0) * x_size * Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! Ay cross entries !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Ay(i,j-1,k-1)
        w = sigy1 * solve_vz / (2.0 * del_z * del_y) - sigz1 * solve_vy / (2.0 * del_z ** 2)
        kk = i + (j - 2) * Nx_nodes + (k - 2) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i,j,k-1)
        w = -sigy0 * solve_vz / (2.0 * del_z * del_y) - sigz1 * solve_vy / (2.0 * del_z ** 2)
        kk = i + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i-1,j-1,k)
        w = -sigx1 * solve_vy / (2.0 * del_x ** 2) + sigy1 * solve_vx / (2.0 * del_x * del_y)
        kk = (i - 1) + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i,j-1,k)
        w = (sigx0 + sigx1) * solve_vy / (2.0 * del_x ** 2) + (sigz0 + sigz1) * solve_vy / (2.0 * del_z ** 2)
        kk = i + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i+1,j-1,k)
        w = -sigx0 * solve_vy / (2.0 * del_x ** 2) - sigy1 * solve_vx / (2.0 * del_x * del_y)
        kk = (i + 1) + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i-1,j,k)
        w = -sigx1 * solve_vy / (2.0 * del_x ** 2) - sigy0 * solve_vx / (2.0 * del_x * del_y)
        kk = (i - 1) + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i,j,k)
        w = (sigx0 + sigx1) * solve_vy / (2.0 * del_x ** 2) + (sigz0 + sigz1) * solve_vy / (2.0 * del_z ** 2)
        kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i+1,j,k)
        w = -sigx0 * solve_vy / (2.0 * del_x ** 2) + sigy0 * solve_vx / (2.0 * del_x * del_y)
        kk = (i + 1) + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i,j-1,k+1)
        w = -sigy1 * solve_vz / (2.0 * del_z * del_y) - sigz0 * solve_vy / (2.0 * del_z ** 2)
        kk = i + (j - 2) * Nx_nodes + (k + 0) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Ay(i,j,k+1)
        w = sigy0 * solve_vz / (2.0 * del_z * del_y) - sigz0 * solve_vy / (2.0 * del_z ** 2)
        kk = i + (j - 1) * Nx_nodes + (k + 0) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! Az cross entries !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Az(i,j-1,k-1)
        w = -sigy1 * solve_vz / (2.0 * del_y ** 2) + sigz1 * solve_vy / (2.0 * del_y * del_z)
        kk = i + (j - 2) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i-1,j,k-1)
        w = -sigx1 * solve_vz / (2.0 * del_x ** 2) + sigz1 * solve_vx / (2.0 * del_x * del_z)
        kk = (i - 1) + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i,j,k-1)
        w = (sigx0 + sigx1) * solve_vz / (2.0 * del_x ** 2) + (sigy0 + sigy1) * solve_vz / (2.0 * del_y ** 2)
        kk = i + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i+1,j,k-1)
        w = -sigx0 * solve_vz / (2.0 * del_x ** 2) - sigz1 * solve_vx / (2.0 * del_x * del_z)
        kk = (i + 1) + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i,j+1,k-1)
        w = -sigy0 * solve_vz / (2.0 * del_y ** 2) - sigz1 * solve_vy / (2.0 * del_y * del_z)
        kk = i + (j + 0) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i,j-1,k)
        w = -sigy1 * solve_vz / (2.0 * del_y ** 2) - sigz0 * solve_vy / (2.0 * del_y * del_z)
        kk = i + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i-1,j,k)
        w = -sigx1 * solve_vz / (2.0 * del_x ** 2) - sigz0 * solve_vx / (2.0 * del_x * del_z)
        kk = (i - 1) + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i,j,k)
        w = (sigx0 + sigx1) * solve_vz / (2.0 * del_x ** 2) + (sigy0 + sigy1) * solve_vz / (2.0 * del_y ** 2)
        kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i+1,j,k)
        w = -sigx0 * solve_vz / (2.0 * del_x ** 2) + sigz0 * solve_vx / (2.0 * del_x * del_z)
        kk = (i + 1) + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !Az(i,j+1,k)
        w = -sigy0 * solve_vz / (2.0 * del_y ** 2) + sigz0 * solve_vy / (2.0 * del_y * del_z)
        kk = i + (j + 0) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1
                                
        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!Lastly V self entries!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        !V(i,j,k-1)
        w = -max(Sigma_z(i,j,k-1), 1.0e-12) / (del_z ** 2)
        kk = i + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !V(i,j-1,k)
        w = -max(Sigma_y(i,j-1,k), 1.0e-12) / (del_y ** 2)
        kk = i + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !V(i-1,j,k)
        w = -max(Sigma_x(i-1,j,k), 1.0e-12) / (del_x ** 2)
        kk = i - 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !V(i,j,k) (Diagonal 7-point Poisson stencil + epsilon floor)
        w = (max(Sigma_x(i,j,k), 1.0e-12) + max(Sigma_x(i-1,j,k), 1.0e-12)) / (del_x ** 2) &
          + (max(Sigma_y(i,j,k), 1.0e-12) + max(Sigma_y(i,j-1,k), 1.0e-12)) / (del_y ** 2) &
          + (max(Sigma_z(i,j,k), 1.0e-12) + max(Sigma_z(i,j,k-1), 1.0e-12)) / (del_z ** 2)
        kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !V(i+1,j,k)
        w = -max(Sigma_x(i,j,k), 1.0e-12) / (del_x ** 2)
        kk = i + 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !V(i,j+1,k)
        w = -max(Sigma_y(i,j,k), 1.0e-12) / (del_y ** 2)
        kk = i + (j + 0) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1

        !V(i,j,k+1)
        w = -max(Sigma_z(i,j,k), 1.0e-12) / (del_z ** 2)
        kk = i + (j - 1) * Nx_nodes + (k + 0) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
        a_mat(nn_counter) = w
        ja(nn_counter) = kk
        nn_counter = nn_counter + 1
        
        !Safety check in case I miscalculated the most possible entries
        if (nn_counter > exact_nn) then
            write(*,*) "Error: a_mat allocation exceeded! Increase multiplier."
            stop
        end if

    end do

    !always ends this way for pardiso
    ia(n_dims + 1) = nn_counter

    !auto check to make sure the ordering hasn't been violated for pardiso
    do rr = 1, n_dims
        do kk = ia(rr), ia(rr+1) - 2
            if (ja(kk) >= ja(kk+1)) then
                write(*,*) "PARDISO Error: Row ", rr, " is unsorted! Col ", ja(kk), " >= ", ja(kk+1)
                stop
            end if
        end do
    end do

    !auto check to make sure there are no negative entries in ja
    if (minval(ja)<0) then
        write(*,*) "PARDISO Error: ja has negative entries"
        stop
    end if

    !Perform Symbolic & Numerical Matrix Factorization and then full solve (Phase 13)
    phase = 13
    call pardiso(pt, maxfct, mnum, mtype, phase, n_dims, a_mat, ia, ja, &
                 i, nrhs, iparm, msglvl, b_vec, x_vec, error)

    !Get A fields back out individually - don't need V
    do k = 1, Nz_nodes
        do j = 1, Ny_nodes
            do i = 1, x_size
                rr = i + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
                Ax(i, j, k) = x_vec(rr)          
            end do
        end do
    end do
    do k = 1, Nz_nodes
        do j = 1, y_size
            do i = 1, Nx_nodes
                rr = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec - Ny_nodes*Nz_nodes
                Ay(i, j, k) = x_vec(rr)          
            end do
        end do
    end do
    do k = 1, z_size
        do j = 1, Ny_nodes
            do i = 1, Nx_nodes
                rr = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes
                Az(i, j, k) = x_vec(rr)          
            end do
        end do
    end do

    ! Compute H field components at 3D cell centers
    do k = 1, z_size
        do j = 1, y_size
            do i = 1, x_size

                ! ==========================================
                ! Hx at cell center
                ! ==========================================
                Hx(i, j, k) = 0.5 * ( &
                    ((Az(i,   j+1, k) - Az(i,   j, k))/del_y - (Ay(i,   j, k+1) - Ay(i,   j, k))/del_z) / Mu_x(i,   j, k) + &
                    ((Az(i+1, j+1, k) - Az(i+1, j, k))/del_y - (Ay(i+1, j, k+1) - Ay(i+1, j, k))/del_z) / Mu_x(i+1, j, k)   &
                )

                ! ==========================================
                ! Hy at cell center 
                ! ==========================================
                Hy(i, j, k) = 0.5 * ( &
                    ((Ax(i, j,   k+1) - Ax(i, j,   k))/del_z - (Az(i+1, j,   k) - Az(i, j,   k))/del_x) / Mu_y(i, j,   k) + &
                    ((Ax(i, j+1, k+1) - Ax(i, j+1, k))/del_z - (Az(i+1, j+1, k) - Az(i, j+1, k))/del_x) / Mu_y(i, j+1, k)   &
                )

                ! ==========================================
                ! Hz at cell center 
                ! ==========================================
                Hz(i, j, k) = 0.5 * ( &
                    ((Ay(i+1, j, k)   - Ay(i, j, k)  )/del_x - (Ax(i, j+1, k)   - Ax(i, j, k)  )/del_y) / Mu_z(i, j, k)   + &
                    ((Ay(i+1, j, k+1) - Ay(i, j, k+1))/del_x - (Ax(i, j+1, k+1) - Ax(i, j, k+1))/del_y) / Mu_z(i, j, k+1) &
                )

            end do
        end do
    end do

    !also save fields, but only keep a 1D slice - for my convenience for now - hard to view 3D fields
    if (mod(counter,tt) == 0) then
        if (TRIM(output_slice)=='x') then
            Hx_out(:,:)=Hx(location_slice,1:y_size,1:z_size)
            Hy_out(:,:)=Hy(location_slice,1:y_size,1:z_size)
            Hz_out(:,:)=Hz(location_slice,1:y_size,1:z_size)
        end if
        if (TRIM(output_slice)=='y') then
            Hx_out(:,:)=Hx(1:x_size,location_slice,1:z_size)
            Hy_out(:,:)=Hy(1:x_size,location_slice,1:z_size)
            Hz_out(:,:)=Hz(1:x_size,location_slice,1:z_size)
        end if
        if (TRIM(output_slice)=='z') then
            Hx_out(:,:)=Hx(1:x_size,1:y_size,location_slice)
            Hy_out(:,:)=Hy(1:x_size,1:y_size,location_slice)
            Hz_out(:,:)=Hz(1:x_size,1:y_size,location_slice)
        end if
        write(3) Hx_out
        write(4) Hy_out
        write(5) Hz_out
    end if

    if (TRIM(save_all_fields)=='yes') then
        write(6) Hx(1:x_size,1:y_size,1:z_size)
        write(7) Hy(1:x_size,1:y_size,1:z_size)
        write(8) Hz(1:x_size,1:y_size,1:z_size)
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Main Transient Step Loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call system_clock(clock_time_end, clock_rate)
    write(*, '("Initialization of steady state completed in ", F6.2, " seconds.")') &
        real(clock_time_end - clock_time_start) / real(clock_rate)
    write(*, '("Starting main time marching algorithm...")')
    call system_clock(clock_time_start)      

    !It's ideal to perform as much computational overhead of pardiso as possible before looping
    !However, the velocity is attached to Ax,Ay,Az so we have to redo phase 2 during each time step
    !Phase 1 can be skipped because ia,ja do not change. But a_mat and b do.
    !If only b was changing, then we wouldn't have to redo step phase 2, just phase 3.
    !So we don't need to update ia,ja anymore.
    !I did add an if statement so that it iwll check if v has changed, if it hasn't at each time step then we can just jump to phase 33 which is fast.

    do counter = 1, time_steps
        
        !at each time step we just need to update a_mat and b_vec, then resolve - phase 23, unless v is not changing

        !There is a small numerical divergence issue from steady to state to first time step if the acceleartion is completely zero
        !I will padd a small acceleration for first time step only to prevent this.

        if (counter==1) then
            if (const_ax*del_t<init_vx*1E-8) then
                temp_acc_x=solve_vx*1E-8/del_t
            end if
            if (const_ay*del_t<init_vy*1E-8) then
                temp_acc_y=solve_vy*1E-8/del_t
            end if
            if (const_az*del_t<init_vz*1E-8) then
                temp_acc_z=solve_vz*1E-8/del_t
            end if
        else
            temp_acc_x=0.0
            temp_acc_y=0.0
            temp_acc_z=0.0
        end if

        solve_vx = trajectory_1(init_vx, const_ax + temp_acc_x, cax_time, counter, del_t)
        solve_vy = trajectory_1(init_vy, const_ay + temp_acc_y, cay_time, counter, del_t)
        solve_vz = trajectory_1(init_vz, const_az + temp_acc_z, caz_time, counter, del_t)

        !for rebuiilding a_mat at each time step
        nn_counter=1

        !!!!!!!!!!!!!!!!!!!!
        !!!!First is Ax!!!!!
        !!!!!!!!!!!!!!!!!!!!
        do rr = 1, n_dims_1vec-Ny_nodes*Nz_nodes
            
            !Ordering here is very important for PARDISO, we need to add them in order down the row
            i = mod(rr - 1, x_size) + 1
            j = mod((rr - 1) / x_size, Ny_nodes) + 1
            k = (rr - 1) / (Ny_nodes * x_size) + 1

            !First check if it's a fixed entry (boundary condition) and if so we know and set A at that position, then move on (cycle)
            if (Fixed_x(i, j, k) == 1.0) then
                x_pos = (i - 1) * del_x
                y_pos = (j - 1) * del_y
                z_pos = (k - 1) * del_z
                b_vec(rr) = 0.5 * By_in * (z_pos) - 0.5 * Bz_in * (y_pos) + &
                1/3.0 * ((gyx*x_pos + gyy*y_pos + gyz*z_pos)*z_pos - (gzx*x_pos + gzy*y_pos + gzz*z_pos)*y_pos)              
                a_mat(nn_counter) = 1.0
                nn_counter = nn_counter + 1
                cycle
            else
                b_vec(rr) = (Sigma_x(i,j,k)/del_t)*Ax(i,j,k)
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!First Ax self entries!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Ax(i,j,k-1)
            w = -1.0 / Mu_y(i,j,k-1) / (del_z ** 2) + (sigma_x(i,j,k) * solve_vz) / (2.0 * del_z)
            kk = i + (j - 1) * x_size + (k - 2) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i,j-1,k)
            w = -1.0 / Mu_z(i,j-1,k) / (del_y ** 2) + (sigma_x(i,j,k) * solve_vy) / (2.0 * del_y)
            kk = i + (j - 2) * x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !GAUGE PENALTY: Ax(i-1,j,k) - Added to enforce -d/dx( 1/mu * dAx/dx )
            w = -1.0 / Mu_x(i-1,j,k) / (del_x ** 2)
            kk = (i - 1) + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i,j,k) (Diagonal spatial + conductivity + gauge penalty)
            w = (1.0 / Mu_z(i,j,k) + 1.0 / Mu_z(i,j-1,k)) / (del_y ** 2) &
            + (1.0 / Mu_y(i,j,k) + 1.0 / Mu_y(i,j,k-1)) / (del_z ** 2) &
            + (1.0 / Mu_x(i,j,k) + 1.0 / Mu_x(i-1,j,k)) / (del_x ** 2) &
            + sigma_x(i,j,k)/del_t
            kk = i + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !GAUGE PENALTY: Ax(i+1,j,k) - Added to enforce -d/dx( 1/mu * dAx/dx )
            w = -1.0 / Mu_x(i,j,k) / (del_x ** 2)
            kk = (i + 1) + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i,j+1,k)
            w = -1.0 / Mu_z(i,j,k) / (del_y ** 2) - (sigma_x(i,j,k) * solve_vy) / (2.0 * del_y)
            kk = i + (j + 0) * x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i,j,k+1)
            w = -1.0 / Mu_y(i,j,k) / (del_z ** 2) - (sigma_x(i,j,k) * solve_vz) / (2.0 * del_z)
            kk = i + (j - 1) * x_size + (k + 0) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!Second Ax has Ay entries
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Ay(i,j-1,k)
            w = 1.0 / Mu_z(i,j-1,k) / (del_x *del_y) - (sigma_x(i,j,k) * solve_vy) / (2.0 * del_x)
            kk = i + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec-Ny_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i+1,j-1,k)
            w = -1.0 / Mu_z(i,j-1,k) / (del_x *del_y) + (sigma_x(i,j,k) * solve_vy) / (2.0 * del_x)
            kk = i + 1 + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec-Ny_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i,j,k)
            w = -1.0 / Mu_z(i,j,k) / (del_x *del_y) - (sigma_x(i,j,k) * solve_vy) / (2.0 * del_x)
            kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec-Ny_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i+1,j,k)
            w = 1.0 / Mu_z(i,j,k) / (del_x *del_y) + (sigma_x(i,j,k) * solve_vy) / (2.0 * del_x)
            kk = i + 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec-Ny_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!Third Ax has Az entries
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Az(i,j,k-1)
            w = 1.0 / Mu_y(i,j,k-1) / (del_x *del_z) - (sigma_x(i,j,k) * solve_vz) / (2.0 * del_x)
            kk = i + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + n_dims_1vec*2-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i+1,j,k-1)
            w = -1.0 / Mu_y(i,j,k-1) / (del_x *del_z) + (sigma_x(i,j,k) * solve_vz) / (2.0 * del_x)
            kk = i + 1 + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + n_dims_1vec*2-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i,j,k)
            w = -1.0 / Mu_y(i,j,k) / (del_x *del_z) - (sigma_x(i,j,k) * solve_vz) / (2.0 * del_x)
            kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i+1,j,k)
            w = 1.0 / Mu_y(i,j,k) / (del_x *del_z) + (sigma_x(i,j,k) * solve_vz) / (2.0 * del_x)
            kk = i + 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!Lastly Ax has V entries
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            !V(i,j,k)
            w = -Sigma_x(i,j,k) / del_x
            kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !V(i+1,j,k)
            w = Sigma_x(i,j,k) / del_x
            kk = i + 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

        end do


        !!!!!!!!!!!!!!!!!!!!
        !!!!Second is Ay!!!!
        !!!!!!!!!!!!!!!!!!!!
        do rr = n_dims_1vec - Ny_nodes*Nz_nodes + 1, n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes

            !Ordering here is very important for PARDISO, we need to add them in order down the row
            rr_rel = rr - (n_dims_1vec - Ny_nodes*Nz_nodes)
            i = mod(rr_rel - 1, Nx_nodes) + 1
            j = mod((rr_rel - 1) / Nx_nodes, y_size) + 1
            k = (rr_rel - 1) / (y_size * Nx_nodes) + 1

            !First check if it's a fixed entry (boundary condition) and if so we know and set A at that position, then move on (cycle)
            if (Fixed_y(i, j, k) == 1.0) then
                x_pos = (i - 1) * del_x
                y_pos = (j - 1) * del_y
                z_pos = (k - 1) * del_z
                b_vec(rr) = 0.5 * Bz_in * (x_pos) - 0.5 * Bx_in * (z_pos) + &
                1/3.0 * ((gzx*x_pos + gzy*y_pos + gzz*z_pos)*x_pos - (gxx*x_pos + gxy*y_pos + gxz*z_pos)*z_pos)                
                a_mat(nn_counter) = 1.0
                nn_counter = nn_counter + 1
                cycle
            else
                b_vec(rr) = (Sigma_y(i,j,k)/del_t)*Ay(i,j,k)
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!First Ay has Ax entries!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Ax(i-1,j,k)
            w = 1.0 / Mu_z(i-1,j,k) / (del_x * del_y) - (Sigma_y(i,j,k) * solve_vx) / (2.0 * del_y)
            kk = i - 1 + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i,j,k)
            w = -1.0 / Mu_z(i,j,k) / (del_x * del_y) - (Sigma_y(i,j,k) * solve_vx) / (2.0 * del_y)
            kk = i + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i-1,j+1,k)
            w = -1.0 / Mu_z(i-1,j,k) / (del_x * del_y) + (Sigma_y(i,j,k) * solve_vx) / (2.0 * del_y)
            kk = i - 1 + (j + 0) * x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i,j+1,k)
            w = 1.0 / Mu_z(i,j,k) / (del_x * del_y) + (Sigma_y(i,j,k) * solve_vx) / (2.0 * del_y)
            kk = i + (j + 0) * x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!Second Ay self entries!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Ay(i,j,k-1)
            w = -1.0 / Mu_x(i,j,k-1) / (del_z ** 2) + (Sigma_y(i,j,k) * solve_vz) / (2.0 * del_z)
            kk = i + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * y_size + n_dims_1vec - Ny_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !GAUGE PENALTY: Ay(i,j-1,k) - Added to enforce -d/dy( 1/mu * dAy/dy )
            w = -1.0 / Mu_y(i,j-1,k) / (del_y ** 2)
            kk = i + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec - Ny_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i-1,j,k)
            w = -1.0 / Mu_z(i-1,j,k) / (del_x ** 2) + (Sigma_y(i,j,k) * solve_vx) / (2.0 * del_x)
            kk = i - 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec - Ny_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i,j,k) (Diagonal spatial + conductivity + gauge penalty)
            w = (1.0 / Mu_z(i,j,k) + 1.0 / Mu_z(i-1,j,k)) / (del_x ** 2) &
                + (1.0 / Mu_x(i,j,k) + 1.0 / Mu_x(i,j,k-1)) / (del_z ** 2) &
                + (1.0 / Mu_y(i,j,k) + 1.0 / Mu_y(i,j-1,k)) / (del_y ** 2) &
                + sigma_y(i,j,k)/del_t
            kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec - Ny_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i+1,j,k)
            w = -1.0 / Mu_z(i,j,k) / (del_x ** 2) - (Sigma_y(i,j,k) * solve_vx) / (2.0 * del_x)
            kk = i + 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec - Ny_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !GAUGE PENALTY: Ay(i,j+1,k) - Added to enforce -d/dy( 1/mu * dAy/dy )
            w = -1.0 / Mu_y(i,j,k) / (del_y ** 2)
            kk = i + (j + 0) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec - Ny_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i,j,k+1)
            w = -1.0 / Mu_x(i,j,k) / (del_z ** 2) - (Sigma_y(i,j,k) * solve_vz) / (2.0 * del_z)
            kk = i + (j - 1) * Nx_nodes + (k + 0) * Nx_nodes * y_size + n_dims_1vec - Ny_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!Third Ay has Az entries!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Az(i,j,k-1)
            w = 1.0 / Mu_x(i,j,k-1) / (del_y * del_z) - (Sigma_y(i,j,k) * solve_vz) / (2.0 * del_y)
            kk = i + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + n_dims_1vec*2-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i,j+1,k-1)
            w = -1.0 / Mu_x(i,j,k-1) / (del_y * del_z) + (Sigma_y(i,j,k) * solve_vz) / (2.0 * del_y)
            kk = i + (j + 0) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + n_dims_1vec*2-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i,j,k)
            w = -1.0 / Mu_x(i,j,k) / (del_y * del_z) - (Sigma_y(i,j,k) * solve_vz) / (2.0 * del_y)
            kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i,j+1,k)
            w = 1.0 / Mu_x(i,j,k) / (del_y * del_z) + (Sigma_y(i,j,k) * solve_vz) / (2.0 * del_y)
            kk = i + (j + 0) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!Lastly Ay has V entries
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            !V(i,j,k)
            w = -Sigma_y(i,j,k) / del_y
            kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !V(i,j+1,k)
            w = Sigma_y(i,j,k) / del_y
            kk = i + (j + 0) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

        end do

        !!!!!!!!!!!!!!!!!!!!
        !!!Third is Az!!!!!!
        !!!!!!!!!!!!!!!!!!!!
        do rr = n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes + 1, n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes

            !Ordering here is very important for PARDISO, we need to add them in order down the row
            rr_rel = rr - (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
            i = mod(rr_rel - 1, Nx_nodes) + 1
            j = mod((rr_rel - 1) / Nx_nodes, Ny_nodes) + 1
            k = (rr_rel - 1) / (Ny_nodes * Nx_nodes) + 1

            !First check if it's a fixed entry (boundary condition) and if so we know and set A at that position, then move on (cycle)
            if (Fixed_z(i, j, k) == 1.0) then
                x_pos = (i - 1) * del_x
                y_pos = (j - 1) * del_y
                z_pos = (k - 1) * del_z
                b_vec(rr) = 0.5 * Bx_in * y_pos - 0.5 * By_in * x_pos + &
                1/3.0 * ((gxx*x_pos + gxy*y_pos + gxz*z_pos)*y_pos - (gyx*x_pos + gyy*y_pos + gyz*z_pos)*x_pos)
                a_mat(nn_counter) = 1.0
                nn_counter = nn_counter + 1
                cycle
            else
                b_vec(rr) = (Sigma_z(i,j,k)/del_t)*Az(i,j,k)
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!First Az has Ax entries!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Ax(i-1,j,k)
            w = 1.0 / Mu_y(i-1,j,k) / (del_x * del_z) - (Sigma_z(i,j,k) * solve_vx) / (2.0 * del_z)
            kk = i - 1 + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i,j,k)
            w = -1.0 / Mu_y(i,j,k) / (del_x * del_z) - (Sigma_z(i,j,k) * solve_vx) / (2.0 * del_z)
            kk = i + (j - 1) *x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i-1,j,k+1)
            w = -1.0 / Mu_y(i-1,j,k) / (del_x * del_z) + (Sigma_z(i,j,k) * solve_vx) / (2.0 * del_z)
            kk = i - 1 + (j - 1) * x_size + (k + 0) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i,j,k+1)
            w = 1.0 / Mu_y(i,j,k) / (del_x * del_z) + (Sigma_z(i,j,k) * solve_vx) / (2.0 * del_z)
            kk = i + (j - 1) * x_size + (k + 0) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!Second Az has Ay entries
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Ay(i,j-1,k)
            w = 1.0 / Mu_x(i,j-1,k) / (del_y * del_z) - (Sigma_z(i,j,k) * solve_vy) / (2.0 * del_z)
            kk = i + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec- Ny_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i,j,k)
            w = -1.0 / Mu_x(i,j,k) / (del_y * del_z) - (Sigma_z(i,j,k) * solve_vy) / (2.0 * del_z)
            kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size +n_dims_1vec- Ny_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i,j-1,k+1)
            w = -1.0 / Mu_x(i,j-1,k) / (del_y * del_z) + (Sigma_z(i,j,k) * solve_vy) / (2.0 * del_z)
            kk = i + (j - 2) * Nx_nodes + (k + 0) * Nx_nodes * y_size + n_dims_1vec- Ny_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i,j,k+1)
            w = 1.0 / Mu_x(i,j,k) / (del_y * del_z) + (Sigma_z(i,j,k) * solve_vy) / (2.0 * del_z)
            kk = i + (j - 1) * Nx_nodes + (k + 0) * Nx_nodes * y_size + n_dims_1vec- Ny_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!Third Az self entries!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            !GAUGE PENALTY: Az(i,j,k-1) - Added to enforce -d/dz( 1/mu * dAz/dz )
            w = -1.0 / Mu_z(i,j,k-1) / (del_z ** 2)
            kk = i + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i,j-1,k)
            w = -1.0 / Mu_x(i,j-1,k) / (del_y ** 2) + (Sigma_z(i,j,k) * solve_vy) / (2.0 * del_y)
            kk = i + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i-1,j,k)
            w = -1.0 / Mu_y(i-1,j,k) / (del_x ** 2) + (Sigma_z(i,j,k) * solve_vx) / (2.0 * del_x)
            kk = i - 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i,j,k) (Diagonal spatial + conductivity + gauge penalty)
            w = (1.0 / Mu_y(i,j,k) + 1.0 / Mu_y(i-1,j,k)) / (del_x ** 2) &
            + (1.0 / Mu_x(i,j,k) + 1.0 / Mu_x(i,j-1,k)) / (del_y ** 2) &
            + (1.0 / Mu_z(i,j,k) + 1.0 / Mu_z(i,j,k-1)) / (del_z ** 2) &
            + sigma_z(i,j,k)/del_t
            kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i+1,j,k)
            w = -1.0 / Mu_y(i,j,k) / (del_x ** 2) - (Sigma_z(i,j,k) * solve_vx) / (2.0 * del_x)
            kk = i + 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i,j+1,k)
            w = -1.0 / Mu_x(i,j,k) / (del_y ** 2) - (Sigma_z(i,j,k) * solve_vy) / (2.0 * del_y)
            kk = i + (j + 0) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !GAUGE PENALTY: Az(i,j,k+1) - Added to enforce -d/dz( 1/mu * dAz/dz )
            w = -1.0 / Mu_z(i,j,k) / (del_z ** 2)
            kk = i + (j - 1) * Nx_nodes + (k + 0) * Nx_nodes * Ny_nodes + n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!Lastly Az has V entries
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            !V(i,j,k)
            w = -Sigma_z(i,j,k) / del_z
            kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !V(i,j,k+1)
            w = Sigma_z(i,j,k) / del_z
            kk = i + (j - 1) * Nx_nodes + (k + 0) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

        end do

        !!!!!!!!!!!!!!!!!!!!
        !!!last is V!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!
        do rr = + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes + 1, n_dims
            
            !Ordering here is very important for PARDISO, we need to add them in order down the row
            rr_rel = rr - (n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes)
            i = mod(rr_rel - 1, Nx_nodes) + 1
            j = mod((rr_rel - 1) / Nx_nodes, Ny_nodes) + 1
            k = (rr_rel - 1) / (Ny_nodes * Nx_nodes) + 1

            !Fix outer nodes to zero for convenience
            if (Fixed_V(i, j, k) == 1.0) then
                x_pos = (i - 1) * del_x
                y_pos = (j - 1) * del_y
                z_pos = (k - 1) * del_z
                b_vec(rr) = - 1.0 * (solve_vy * Bz_in - solve_vz * By_in) * x_pos &
                - 1.0 * (solve_vz * Bx_in - solve_vx * Bz_in) * y_pos &
                - 1.0 * (solve_vx * By_in - solve_vy * Bx_in) * z_pos &
                - 0.5 * (solve_vy * (gzx*x_pos + gzy*y_pos + gzz*z_pos) - solve_vz * (gyx*x_pos + gyy*y_pos + gyz*z_pos)) * x_pos &
                - 0.5 * (solve_vz * (gxx*x_pos + gxy*y_pos + gxz*z_pos) - solve_vx * (gzx*x_pos + gzy*y_pos + gzz*z_pos)) * y_pos &
                - 0.5 * (solve_vx * (gyx*x_pos + gyy*y_pos + gyz*z_pos) - solve_vy * (gxx*x_pos + gxy*y_pos + gxz*z_pos)) * z_pos
                a_mat(nn_counter) = 1.0
                nn_counter = nn_counter + 1
                cycle
            else
            !Divergence of (sigma * A / del_t) from the previous time step
                b_vec(rr) = &
                    ( Sigma_x(i,j,k) * Ax(i,j,k) - Sigma_x(i-1,j,k) * Ax(i-1,j,k) ) / (del_t * del_x) + &
                    ( Sigma_y(i,j,k) * Ay(i,j,k) - Sigma_y(i,j-1,k) * Ay(i,j-1,k) ) / (del_t * del_y) + &
                    ( Sigma_z(i,j,k) * Az(i,j,k) - Sigma_z(i,j,k-1) * Az(i,j,k-1) ) / (del_t * del_z)
            end if


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!! div( sigma * (v x curl A) ) term added to the V equation !!!
            !!! (fully merged across x/y/z directions, sorted by kk)      !!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            sigx0 = Sigma_x(i,j,k)
            sigx1 = Sigma_x(i-1,j,k)
            sigy0 = Sigma_y(i,j,k)
            sigy1 = Sigma_y(i,j-1,k)
            sigz0 = Sigma_z(i,j,k)
            sigz1 = Sigma_z(i,j,k-1)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!! Ax cross entries !!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Ax(i-1,j,k-1)
            w = sigx1 * solve_vz / (2.0 * del_z * del_x) - sigz1 * solve_vx / (2.0 * del_z ** 2)
            kk = (i - 1) + (j - 1) * x_size + (k - 2) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i,j,k-1)
            w = -sigx0 * solve_vz / (2.0 * del_z * del_x) - sigz1 * solve_vx / (2.0 * del_z ** 2)
            kk = i + (j - 1) * x_size + (k - 2) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i-1,j-1,k)
            w = sigx1 * solve_vy / (2.0 * del_y * del_x) - sigy1 * solve_vx / (2.0 * del_y ** 2)
            kk = (i - 1) + (j - 2) * x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i,j-1,k)
            w = -sigx0 * solve_vy / (2.0 * del_y * del_x) - sigy1 * solve_vx / (2.0 * del_y ** 2)
            kk = i + (j - 2) * x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            ja(nn_counter) = kk
            nn_counter = nn_counter + 1

            !Ax(i-1,j,k)
            w = (sigy0 + sigy1) * solve_vx / (2.0 * del_y ** 2) + (sigz0 + sigz1) * solve_vx / (2.0 * del_z ** 2) &
            - sigx1 / (del_x *del_t)
            kk = (i - 1) + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i,j,k)
            w = (sigy0 + sigy1) * solve_vx / (2.0 * del_y ** 2) + (sigz0 + sigz1) * solve_vx / (2.0 * del_z ** 2) &
            + sigx0 / (del_x *del_t)
            kk = i + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i-1,j+1,k)
            w = -sigx1 * solve_vy / (2.0 * del_y * del_x) - sigy0 * solve_vx / (2.0 * del_y ** 2)
            kk = (i - 1) + (j + 0) * x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i,j+1,k)
            w = sigx0 * solve_vy / (2.0 * del_y * del_x) - sigy0 * solve_vx / (2.0 * del_y ** 2)
            kk = i + (j + 0) * x_size + (k - 1) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i-1,j,k+1)
            w = -sigx1 * solve_vz / (2.0 * del_z * del_x) - sigz0 * solve_vx / (2.0 * del_z ** 2)
            kk = (i - 1) + (j - 1) * x_size + (k + 0) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ax(i,j,k+1)
            w = sigx0 * solve_vz / (2.0 * del_z * del_x) - sigz0 * solve_vx / (2.0 * del_z ** 2)
            kk = i + (j - 1) * x_size + (k + 0) * x_size * Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!! Ay cross entries !!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Ay(i,j-1,k-1)
            w = sigy1 * solve_vz / (2.0 * del_z * del_y) - sigz1 * solve_vy / (2.0 * del_z ** 2)
            kk = i + (j - 2) * Nx_nodes + (k - 2) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i,j,k-1)
            w = -sigy0 * solve_vz / (2.0 * del_z * del_y) - sigz1 * solve_vy / (2.0 * del_z ** 2)
            kk = i + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i-1,j-1,k)
            w = -sigx1 * solve_vy / (2.0 * del_x ** 2) + sigy1 * solve_vx / (2.0 * del_x * del_y)
            kk = (i - 1) + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i,j-1,k)
            w = (sigx0 + sigx1) * solve_vy / (2.0 * del_x ** 2) + (sigz0 + sigz1) * solve_vy / (2.0 * del_z ** 2) &
            - sigy1 / (del_y * del_t)
            kk = i + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i+1,j-1,k)
            w = -sigx0 * solve_vy / (2.0 * del_x ** 2) - sigy1 * solve_vx / (2.0 * del_x * del_y)
            kk = (i + 1) + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i-1,j,k)
            w = -sigx1 * solve_vy / (2.0 * del_x ** 2) - sigy0 * solve_vx / (2.0 * del_x * del_y)
            kk = (i - 1) + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i,j,k)
            w = (sigx0 + sigx1) * solve_vy / (2.0 * del_x ** 2) + (sigz0 + sigz1) * solve_vy / (2.0 * del_z ** 2) &
            + sigy0 / (del_y * del_t)
            kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i+1,j,k)
            w = -sigx0 * solve_vy / (2.0 * del_x ** 2) + sigy0 * solve_vx / (2.0 * del_x * del_y)
            kk = (i + 1) + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i,j-1,k+1)
            w = -sigy1 * solve_vz / (2.0 * del_z * del_y) - sigz0 * solve_vy / (2.0 * del_z ** 2)
            kk = i + (j - 2) * Nx_nodes + (k + 0) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Ay(i,j,k+1)
            w = sigy0 * solve_vz / (2.0 * del_z * del_y) - sigz0 * solve_vy / (2.0 * del_z ** 2)
            kk = i + (j - 1) * Nx_nodes + (k + 0) * Nx_nodes * y_size + (n_dims_1vec - Ny_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!! Az cross entries !!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Az(i,j-1,k-1)
            w = -sigy1 * solve_vz / (2.0 * del_y ** 2) + sigz1 * solve_vy / (2.0 * del_y * del_z)
            kk = i + (j - 2) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i-1,j,k-1)
            w = -sigx1 * solve_vz / (2.0 * del_x ** 2) + sigz1 * solve_vx / (2.0 * del_x * del_z)
            kk = (i - 1) + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i,j,k-1)
            w = (sigx0 + sigx1) * solve_vz / (2.0 * del_x ** 2) + (sigy0 + sigy1) * solve_vz / (2.0 * del_y ** 2) &
            - sigz1 / (del_z * del_t)
            kk = i + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i+1,j,k-1)
            w = -sigx0 * solve_vz / (2.0 * del_x ** 2) - sigz1 * solve_vx / (2.0 * del_x * del_z)
            kk = (i + 1) + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i,j+1,k-1)
            w = -sigy0 * solve_vz / (2.0 * del_y ** 2) - sigz1 * solve_vy / (2.0 * del_y * del_z)
            kk = i + (j + 0) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i,j-1,k)
            w = -sigy1 * solve_vz / (2.0 * del_y ** 2) - sigz0 * solve_vy / (2.0 * del_y * del_z)
            kk = i + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i-1,j,k)
            w = -sigx1 * solve_vz / (2.0 * del_x ** 2) - sigz0 * solve_vx / (2.0 * del_x * del_z)
            kk = (i - 1) + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i,j,k)
            w = (sigx0 + sigx1) * solve_vz / (2.0 * del_x ** 2) + (sigy0 + sigy1) * solve_vz / (2.0 * del_y ** 2) &
            + sigz0 / (del_z * del_t)
            kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i+1,j,k)
            w = -sigx0 * solve_vz / (2.0 * del_x ** 2) + sigz0 * solve_vx / (2.0 * del_x * del_z)
            kk = (i + 1) + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !Az(i,j+1,k)
            w = -sigy0 * solve_vz / (2.0 * del_y ** 2) + sigz0 * solve_vy / (2.0 * del_y * del_z)
            kk = i + (j + 0) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + (n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes)
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1
                                                        
            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!Lastly V self entries!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            !V(i,j,k-1)
            w = -max(Sigma_z(i,j,k-1), 1.0e-12) / (del_z ** 2)
            kk = i + (j - 1) * Nx_nodes + (k - 2) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !V(i,j-1,k)
            w = -max(Sigma_y(i,j-1,k), 1.0e-12) / (del_y ** 2)
            kk = i + (j - 2) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !V(i-1,j,k)
            w = -max(Sigma_x(i-1,j,k), 1.0e-12) / (del_x ** 2)
            kk = i - 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !V(i,j,k) (Diagonal 7-point Poisson stencil + epsilon floor)
            w = (max(Sigma_x(i,j,k), 1.0e-12) + max(Sigma_x(i-1,j,k), 1.0e-12)) / (del_x ** 2) &
            + (max(Sigma_y(i,j,k), 1.0e-12) + max(Sigma_y(i,j-1,k), 1.0e-12)) / (del_y ** 2) &
            + (max(Sigma_z(i,j,k), 1.0e-12) + max(Sigma_z(i,j,k-1), 1.0e-12)) / (del_z ** 2)
            kk = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !V(i+1,j,k)
            w = -max(Sigma_x(i,j,k), 1.0e-12) / (del_x ** 2)
            kk = i + 1 + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !V(i,j+1,k)
            w = -max(Sigma_y(i,j,k), 1.0e-12) / (del_y ** 2)
            kk = i + (j + 0) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1

            !V(i,j,k+1)
            w = -max(Sigma_z(i,j,k), 1.0e-12) / (del_z ** 2)
            kk = i + (j - 1) * Nx_nodes + (k + 0) * Nx_nodes * Ny_nodes + n_dims_1vec*3-Ny_nodes*Nz_nodes-Nx_nodes*Nz_nodes-Nx_nodes*Ny_nodes
            a_mat(nn_counter) = w
            nn_counter = nn_counter + 1
        
        end do

        if ((o_solve_vx == solve_vx .and. o_solve_vy == solve_vy .and. o_solve_vz == solve_vz) .and. (counter > 1)) then
            !Direct back-substitution step solve phase (Phase 33)
            phase = 33
            call pardiso(pt, maxfct, mnum, mtype, phase, n_dims, a_mat, ia, ja, &
                        i, nrhs, iparm, msglvl, b_vec, x_vec, error)
        else
            !Factor and Direct back-substitution step solve phase (Phase 23)
            phase = 23
            call pardiso(pt, maxfct, mnum, mtype, phase, n_dims, a_mat, ia, ja, &
                        i, nrhs, iparm, msglvl, b_vec, x_vec, error)
        end if

        !udpate o_solvers for next round
        o_solve_vx = solve_vx
        o_solve_vy = solve_vy
        o_solve_vz = solve_vz

        !Get A fields back out individually at each time step - don't need V
        do k = 1, Nz_nodes
            do j = 1, Ny_nodes
                do i = 1, x_size
                    rr = i + (j - 1) * x_size + (k - 1) * x_size * Ny_nodes
                    Ax(i, j, k) = x_vec(rr)          
                end do
            end do
        end do
        do k = 1, Nz_nodes
            do j = 1, y_size
                do i = 1, Nx_nodes
                    rr = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * y_size + n_dims_1vec - Ny_nodes*Nz_nodes
                    Ay(i, j, k) = x_vec(rr)          
                end do
            end do
        end do
        do k = 1, z_size
            do j = 1, Ny_nodes
                do i = 1, Nx_nodes
                    rr = i + (j - 1) * Nx_nodes + (k - 1) * Nx_nodes * Ny_nodes + n_dims_1vec*2 - Ny_nodes*Nz_nodes - Nx_nodes*Nz_nodes
                    Az(i, j, k) = x_vec(rr)          
                end do
            end do
        end do

        ! Compute H field components at 3D cell centers
        do k = 1, z_size
            do j = 1, y_size
                do i = 1, x_size

                    ! ==========================================
                    ! Hx at cell center
                    ! ==========================================
                    Hx(i, j, k) = 0.5 * ( &
                        ((Az(i,   j+1, k) - Az(i,   j, k))/del_y - (Ay(i,   j, k+1) - Ay(i,   j, k))/del_z) / Mu_x(i,   j, k) + &
                        ((Az(i+1, j+1, k) - Az(i+1, j, k))/del_y - (Ay(i+1, j, k+1) - Ay(i+1, j, k))/del_z) / Mu_x(i+1, j, k)   &
                    )

                    ! ==========================================
                    ! Hy at cell center 
                    ! ==========================================
                    Hy(i, j, k) = 0.5 * ( &
                        ((Ax(i, j,   k+1) - Ax(i, j,   k))/del_z - (Az(i+1, j,   k) - Az(i, j,   k))/del_x) / Mu_y(i, j,   k) + &
                        ((Ax(i, j+1, k+1) - Ax(i, j+1, k))/del_z - (Az(i+1, j+1, k) - Az(i, j+1, k))/del_x) / Mu_y(i, j+1, k)   &
                    )

                    ! ==========================================
                    ! Hz at cell center 
                    ! ==========================================
                    Hz(i, j, k) = 0.5 * ( &
                        ((Ay(i+1, j, k)   - Ay(i, j, k)  )/del_x - (Ax(i, j+1, k)   - Ax(i, j, k)  )/del_y) / Mu_z(i, j, k)   + &
                        ((Ay(i+1, j, k+1) - Ay(i, j, k+1))/del_x - (Ax(i, j+1, k+1) - Ax(i, j, k+1))/del_y) / Mu_z(i, j, k+1) &
                    )

                end do
            end do
        end do

        !writes out what time step we are out so we can track it
        if (mod(counter,10) == 0) then
            write(*, '(I0, " of ", I0, " time steps")') counter, time_steps
        end if

        !also save fields at user designated time steps, but only keep a 1D slice - for my convenience for now - hard to view 3D fields
        if (mod(counter,tt) == 0) then
            if (TRIM(output_slice)=='x') then
                Hx_out(:,:)=Hx(location_slice,1:y_size,1:z_size)
                Hy_out(:,:)=Hy(location_slice,1:y_size,1:z_size)
                Hz_out(:,:)=Hz(location_slice,1:y_size,1:z_size)
            end if
            if (TRIM(output_slice)=='y') then
                Hx_out(:,:)=Hx(1:x_size,location_slice,1:z_size)
                Hy_out(:,:)=Hy(1:x_size,location_slice,1:z_size)
                Hz_out(:,:)=Hz(1:x_size,location_slice,1:z_size)
            end if
            if (TRIM(output_slice)=='z') then
                Hx_out(:,:)=Hx(1:x_size,1:y_size,location_slice)
                Hy_out(:,:)=Hy(1:x_size,1:y_size,location_slice)
                Hz_out(:,:)=Hz(1:x_size,1:y_size,location_slice)
            end if
            write(3) Hx_out
            write(4) Hy_out
            write(5) Hz_out
        end if

        if (TRIM(save_all_fields)=='yes') then
            write(6) Hx(1:x_size,1:y_size,1:z_size)
            write(7) Hy(1:x_size,1:y_size,1:z_size)
            write(8) Hz(1:x_size,1:y_size,1:z_size)
        end if

    end do

    !close the binary read outs since done with time looping
    close(3)
    close(4)
    close(5)

    if (TRIM(save_all_fields)=='yes') then
        close(6)
        close(7)
        close(8)
    end if
    
    !Completely clear PARDISO internal structural handles after time looping
    phase = -1
    call pardiso(pt, maxfct, mnum, mtype, phase, n_dims, a_mat, ia, ja, &
                 i, nrhs, iparm, msglvl, b_vec, x_vec, error)

    call system_clock(clock_time_end, clock_rate)
    write(*, '("Main time marching completed in ", F6.2, " seconds.")') &
        real(clock_time_end - clock_time_start) / real(clock_rate)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! functions & utilities !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    contains

        Function trajectory_1(v_fx,a_fx,t_fx,current_timestep_fx,delt_fx)
            implicit none
            real(dp) trajectory_1
            real(dp) :: v_fx, a_fx, delt_fx
            integer :: current_timestep_fx, t_fx
            
            !v_fx is initial velocity
            !a_fx if the acceleration
            !t_fx is the number of time steps to apply acceleration, after which v is constant based on last acceleration
            !delt_fx is the time step delta
            !current_timestep_fx is the current time step integer

            if (current_timestep_fx<=t_fx) then
                trajectory_1=v_fx+a_fx*current_timestep_fx*delt_fx
            end if
            if (current_timestep_fx>t_fx) then
                trajectory_1=v_fx+a_fx*t_fx*delt_fx
            end if

        End Function trajectory_1

End Program Quasimagnetostaticsolver3D