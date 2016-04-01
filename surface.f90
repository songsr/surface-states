
 
  !==================================================================!
  program surface_stat
  !==================================================================!




  implicit none


  integer, parameter          :: dp = selected_real_kind(15,300)
  
  real(kind=dp), parameter    :: pi=3.141592653589793238462643383279_dp
  real(kind=dp), parameter    :: twopi = 2*pi
  complex(kind=dp), parameter :: cmplx_i = (0.0_dp,1.0_dp)
  complex(kind=dp), parameter :: cmplx_0 = (0.0_dp,0.0_dp)
  complex(kind=dp), parameter :: cmplx_1 = (1.0_dp,0.0_dp)

  real(kind=dp), parameter :: bohr = 0.5291772108_dp

  real(kind=dp), parameter    :: eps2  = 1.0e-2_dp
  real(kind=dp), parameter    :: eps5  = 1.0e-5_dp
  real(kind=dp), parameter    :: eps6  = 1.0e-6_dp
  real(kind=dp), parameter    :: eps7  = 1.0e-7_dp
  real(kind=dp), parameter    :: eps8  = 1.0e-8_dp
  real(kind=dp), parameter    :: eps10 = 1.0e-10_dp
  
 
  
  integer:: num_wann
  integer:: mp_grid(3)
  real(kind=dp) :: dist_cutoff
  real(kind=dp) :: hr_cutoff
  real(kind=dp) :: fermi_energy
  real(kind=dp) :: tran_win_min
  real(kind=dp) :: tran_win_max
  real(kind=dp) :: tran_energy_step
  real(kind=dp) :: real_lattice(3,3)



 
  ! Hamiltonian matrix in WF representation
  !
  complex(kind=dp), allocatable :: ham_r(:,:,:)
  !
  ! irvec(i,irpt)     The irpt-th Wigner-Seitz grid point has components
  !                   irvec(1:3,irpt) in the basis of the lattice vectors
  !
  integer,  allocatable :: irvec(:,:)
  !
  ! ndegen(irpt)      Weight of the irpt-th point is 1/ndegen(irpt)
  !
  integer,   allocatable :: ndegen(:)
  !
  ! nrpts             number of Wigner-Seitz grid points
  !
  integer:: nrpts
  !
  ! translated Wannier centres
  !
  real(kind=dp), allocatable :: wannier_centres_translated(:,:) 

    
  
! small complex number 
  complex(kind=dp), parameter :: eta=(0.d0,0.0005d0)

! nterx  = # of maximum iteration to calculate transfer matrix
  integer, parameter :: nterx=50
  ! cartesian axis to which real_lattice(:,one_dim_vec) is parallel
  integer :: one_dim_vec
  integer :: one_dim_dir
  integer :: nrpts_one_dim
  ! num_pl : number of unit cell in a principal layer
  integer :: num_pl
  
      integer ::   tran_num_bb
      
  ! coord : coord(1) defines the conduction direction according to 1=x,2=y,3=z, 
  ! coord(2),coord(3) define the other directions during sorting routines
!  integer,dimension(3) :: coord
  ! index of sorted WF centres to unsorted
 ! integer,allocatable :: tran_sorted_idx(:)

!  real(kind=dp), allocatable :: hr_one_dim(:,:,:)
!  real(kind=dp), allocatable :: hB0(:,:)
!  real(kind=dp), allocatable :: hB1(:,:)


  complex(kind=dp), allocatable :: hr_one_dim(:,:,:)
  complex(kind=dp), allocatable :: hB0(:,:)
  complex(kind=dp), allocatable :: hB1(:,:)

!   integer stdout
    integer ierr
    integer file_unit




    real(kind=dp) k    
    real(kind=dp) e_scan
    integer n,ikp,n_e
    integer i,j

     real(dp), allocatable :: omega(:)
     
     real(dp), allocatable :: dos_l(:,:)
     real(dp), allocatable :: dos_r(:,:)    
          real(dp), allocatable :: dos(:,:)

     integer knv2
     character*12 :: fname='input.dat'
     logical ::  exists
     
     
     knv2=299
     
     mp_grid(1)=12
     mp_grid(2)=12
     mp_grid(3)=12
     
     one_dim_vec=2
     one_dim_dir=2
     tran_energy_step=0.01
     tran_win_max=-0.5
     tran_win_min=-11.0
     
     
     do i=1,3
         do j=1,3
             real_lattice(i,j)=0
         enddo
     enddo
     
     do i=1,3
         do j=1,3
             wannier_centres_translated(j,i)=0
         enddo
     enddo     
     
     
     
      inquire(file=fname,exist=exists)
     if (exists)then
        write(*,*) 'read some paramters from input.dat'
        open(unit=1001,file=fname,status='old')
     else
         write(*,*)'file' ,fname, 'dos not exist'
        stop
     endif
 
      read(1001,*) num_wann
      read(1001,*)knv2
      read(1001,*)mp_grid(:)     
       read(1001,*)one_dim_vec     
      read(1001,*)one_dim_dir     
      read(1001,*)tran_energy_step     
      read(1001,*)tran_win_max 
      read(1001,*)tran_win_min 
 
     do i=1,3
         read(1001,*) real_lattice(:,i)
     enddo
     
     do i=1,num_wann
              read(1001,*)  wannier_centres_translated(:,i)
     enddo     
     
     

     
    write(*,'(/1x,a)') '*---------------------------------------------------------------------------*'
    write(*,'(1x,a)') '|                              Surface Dos                                    |'
    write(*,'(1x,a)') '*---------------------------------------------------------------------------*'
    write(*,*)


      write(*,'(/1x,a/)') 'Calculation of Surface Density of States'

     
          call hamiltonian_setup()
          call hamiltonian_read_hr()
    !   loop over the energies

    n_e = floor((tran_win_max-tran_win_min)/tran_energy_step)+1

     allocate( dos(knv2, n_e))
     
     allocate( omega(n_e))
     allocate( dos_l(knv2, n_e))
     allocate( dos_r(knv2, n_e))
!    write(*,'(/1x,a)',advance='no') 'Calculating quantum&
!        & conductance and density of states...'

        open (unit=102, file='dos.dat_l')

   do ikp=1,300,1 !N_kpoint
      k=-0.5 + (ikp-1)*1.0/real(knv2-1,dp)
     do n=1,n_e
       e_scan = tran_win_min + real(n-1,dp)*tran_energy_step

          call reduce_hr(k)
          call cut_hr_one_dim()
          call get_ht()

          call surface_dos(e_scan,dos(ikp,n))
          write(102, '(3f16.8)') k, e_scan, dos(ikp, n)   
          deallocatable(hr_one_dim, hB0, hB1)
      enddo
      write(102,*)
    end do
    
    close(102)
    
    
    

 
 
 contains
   !============================================!
  subroutine hamiltonian_setup()
    !============================================!

    implicit none
    character(len=10) header
    file_unit=11
    
    open(file_unit,file='wannier90_hr.dat',form='formatted',status='unknown',err=101)


    read(file_unit,*) header ! Date and time
    read(file_unit,*) num_wann
    read(file_unit,*) nrpts
    close(file_unit)

 
    allocate(irvec(3,nrpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating irvec in hamiltonian_setup')
    irvec=0
    !
    allocate(ndegen(nrpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ndegen in hamiltonian_setup')
    ndegen=0
    !
    allocate(ham_r(num_wann,num_wann,nrpts),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ham_r in hamiltonian_setup')
    ham_r=cmplx_0

    !
    allocate(wannier_centres_translated(3,num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error allocating wannier_centres_translated in hamiltonian_setup')
    wannier_centres_translated=0.0_dp
    
    return
101 call io_error('Error: Problem opening input file wannier90_hr.dat')
  end subroutine hamiltonian_setup
  

   !============================================!
  subroutine hamiltonian_read_hr()
    !============================================!
    !  Read the Hamiltonian in the WF basis     !
    !============================================!


    integer            :: i,j,loop_rpt,file_unit,i1,j1
    character (len=33) :: header


    ! read the  whole matrix with all the indices 
 
    file_unit=11
    open(file_unit,file='wannier90_hr.dat',form='formatted',status='unknown',err=101)


    read(file_unit,*) header ! Date and time
    read(file_unit,*) num_wann
    read(file_unit,*) nrpts
 
    read(file_unit,'(15I5)') (ndegen(i),i=1,nrpts)
    do loop_rpt=1,nrpts
       do i=1,num_wann
          do j=1,num_wann
             read( file_unit,'(5I5,2F12.6)') irvec(:,loop_rpt), j1, i1, ham_r(j,i,loop_rpt)
          end do
       end do
    end do

    close(file_unit)

    return

101 call io_error('Error: Problem opening input file wannier90_hr.dat')

  end subroutine hamiltonian_read_hr 
  
  
 
    !==================================================================!
       subroutine io_error ( error_msg )
    !==================================================================!
    !                                                                  !
    ! Aborts giving error message                                      !
    !                                                                  !
    !===================================================================  

         implicit none
         character(len=*), intent(in) :: error_msg

         write(*,*)  'Exiting.......' 
         write(*, '(1x,a)') trim(error_msg)
         
 !        close(*)
         
         stop "wannier90 error: examine the output/error file for details" 
         
       end subroutine io_error 
       
       
       
 
    !==================================================================!
  subroutine reduce_hr(kx)
    !==================================================================!
    !
    ! reduce ham_r from 3-d to 1-d
    !

    implicit none

! wave vector k times lattice vector R  
     real(Dp) :: kdotr
  complex(kind=dp)   ratio 
! input wave vector k's cooridinates
     real(Dp),intent(in) :: kx  !1 dimension edge 
    integer :: ierr
    integer :: irvec_max, irvec_tmp(3), two_dim_vec(2)
    integer :: i, j
    integer :: i1, i2, i3, n1, nrpts_tmp, loop_rpt
 

    one_dim_vec=2  ! y direction,semi-infinite direction
       
    j=0
    do i=1,3
       if ( i .ne. one_dim_vec ) then
          j = j +1
          two_dim_vec(j)=i
       end if
    end do

    ! starting H matrix should include all W-S supercell where 
    ! the center of the cell spans the full space of the home cell
    ! adding one more buffer layer when mp_grid(one_dim_vec) is an odd number

    !irvec_max = (mp_grid(one_dim_vec)+1)/2
    irvec_tmp = maxval(irvec,DIM=2)+1
    irvec_max = irvec_tmp(one_dim_vec)
    nrpts_one_dim = 2*irvec_max+1
 
    allocate(hr_one_dim(num_wann,num_wann,-irvec_max:irvec_max),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating hr_one_dim in tran_reduce_hr')  
    hr_one_dim = 0.0_dp
       
    ! check imaginary part
!    write(*,'(1x,a,F12.6)') 'Maximum imaginary part of the real-space Hamiltonian: ',maxval(abs(aimag(ham_r)))

    ! select a subset of ham_r, where irvec is 0 along the Z lattice vector

    nrpts_tmp = 0
loop_n1: do n1 = -irvec_max, irvec_max
       do loop_rpt = 1,nrpts
          i1 = mod(n1 - irvec(one_dim_vec,loop_rpt),mp_grid(one_dim_vec))
          i2 = irvec(two_dim_vec(1),loop_rpt)
          i3 = irvec(two_dim_vec(2),loop_rpt)
          if (i1.eq.0  .and. i3.eq.0 ) then
             
             kdotr=kx*i2
             ratio=cos(twopi*kdotr)+cmplx_i*sin(twopi*kdotr)             
             hr_one_dim(:,:,n1) = hr_one_dim(:,:,n1)+ ham_r(:,:,loop_rpt)*ratio/ndegen(loop_rpt)
            
          end if
       end do
    end do loop_n1   


    return

  end subroutine reduce_hr
  
  
  
   !==================================================================!
  subroutine cut_hr_one_dim()
    !==================================================================!
    !

    implicit none
    !
    integer :: irvec_max
    integer :: i, j, n1
    real(kind=dp) :: hr_max
    real(kind=dp) :: dist
    real(kind=dp) :: dist_vec(3)
    real(kind=dp) :: dist_ij_vec(3)
    real(kind=dp) :: shift_vec(3,-nrpts_one_dim/2:nrpts_one_dim/2)
    real(kind=dp) :: hr_tmp(num_wann,num_wann)
    character(len=20) :: length_unit='ang'
 
 
    irvec_max = nrpts_one_dim/2
    ! maximum possible dist_cutoff
    dist = real(mp_grid(one_dim_vec),dp)*abs(real_lattice(one_dim_dir,one_dim_vec))/2.0_dp

    if ( dist_cutoff .gt. dist ) then
       write(*,'(1x,a,1x,F10.5,1x,a)') 'dist_cutoff',dist_cutoff,trim(length_unit),'is too large'
       dist_cutoff = dist
       write(*,'(4x,a,1x,F10.5,1x,a)') 'reset to',dist_cutoff,trim(length_unit)
    end if

    do n1 = -irvec_max, irvec_max
       shift_vec(:,n1) = real(n1,dp)*(real_lattice(:,one_dim_vec))
!           write(*,'(a,3f10.6)') 'shift_vec', shift_vec(:,n1)
    end do

    ! apply dist_cutoff first

       do i=1,num_wann
          do j=1,num_wann
             dist_ij_vec(one_dim_dir)=wannier_centres_translated(one_dim_dir,i)-wannier_centres_translated(one_dim_dir,j)
             do n1 = -irvec_max, irvec_max
                dist_vec(one_dim_dir) = dist_ij_vec(one_dim_dir)+ shift_vec(one_dim_dir,n1)
                dist = abs(dist_vec(one_dim_dir))
                if ( dist .gt. dist_cutoff ) hr_one_dim(j,i,n1)=0.0_dp
             end do
          end do
        end do


    ! output maximum to check a decay of H as a function of lattice vector R
!    write(*,'(/1x,a78)') repeat('-',78)
!    write(*,'(1x,4x,a)') &
!                'Maximum real part of the real-space Hamiltonian at each lattice point'
!    write(*,'(1x,8x,a62)') repeat('-',62)
!    write(*,'(1x,11x,a,11x,a)') 'Lattice point R', 'Max |H_ij(R)|'
    ! calculate number of units inside a principal layer
    num_pl = 0
    do n1=-irvec_max,irvec_max
       hr_tmp(:,:)=abs(hr_one_dim(:,:,n1))       
       hr_max = maxval(hr_tmp) 
       if ( hr_max .gt. hr_cutoff ) then
          if (abs(n1) .gt. num_pl) num_pl = abs(n1)
       else
          hr_one_dim(:,:,n1)=0.0_dp
       end if
       write(*,'(1x,9x,5x,I5,5x,12x,F12.6)') n1, hr_max 
    end do
!    write(*,'(1x,8x,a62)') repeat('-',62)

!        write(*,'(/1x,a,I6)') 'Number of unit cells inside the principal layer:',num_pl 
!        write(*,'(1x,a,I6)')  'Number of Wannier Functions inside the principal layer:',num_pl*num_wann 

    ! apply hr_cutoff to each element inside the principal layer
    do n1 = -num_pl , num_pl
       do i=1,num_wann
          do j=1,num_wann
             if ( abs(hr_one_dim(j,i,n1)) .lt. hr_cutoff ) hr_one_dim(j,i,n1)=0.0_dp
          end do
       end do
    end do



    return
  
  end subroutine cut_hr_one_dim
  
  
  
  
  
  
  !==================================================================!
    subroutine get_ht()
    !==================================================================!
    !  construct h00 and h01
    !==================================================================!
    !

    !
    implicit none
    !
    integer :: ierr, file_unit
    integer :: i, j, n1, im, jm
    character(len=9)   :: cdate, ctime

    !

    !
    tran_num_bb=num_pl*num_wann
    !
    allocate(hB0(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating hB0 in tran_get_ht')
    allocate(hB1(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating hB1 in tran_get_ht')
    !
    hB0 = 0.0_dp 
    hB1 = 0.0_dp 
    !
    ! h00
    do j=0, num_pl-1
       do i=0, num_pl-1
          n1 = i-j
          im=i*num_wann
          jm=j*num_wann
          hB0(jm+1:jm+num_wann,im+1:im+num_wann)=hr_one_dim(:,:,n1)
       end do
    end do
  
    ! h01
    do j=1, num_pl
       do i=0, j-1
          n1 = i-(j-1)+num_pl
          im=i*num_wann
          jm=(j-1)*num_wann
          hB1(jm+1:jm+num_wann,im+1:im+num_wann)=hr_one_dim(:,:,n1)
       end do
    end do

    ! shift by fermi_energy
    do i=1,tran_num_bb
       hB0(i,i)=hB0(i,i)-fermi_energy
    end do


  
!       file_unit = 13
!       open(file_unit,file='wannier90_htB.dat',status='unknown',form='formatted',action='write')

!       call io_date(cdate,ctime)
!       write(file_unit,*) 'written on '//cdate//' at '//ctime ! Date and time
!       write(file_unit,'(I6)') tran_num_bb
!       write(file_unit,'(6F12.6)') ((hB0(j,i),j=1,tran_num_bb),i=1,tran_num_bb)
!       write(file_unit,'(I6)') tran_num_bb
!       write(file_unit,'(6F12.6)') ((hB1(j,i),j=1,tran_num_bb),i=1,tran_num_bb)
     
!       close(file_unit)
 
 

    return
  
  end subroutine get_ht
  
  
  
  
   !==================================================================!
  subroutine surface_dos(e_scan,dos)
    !==================================================================!

    implicit none

    integer :: qc_unit, dos_unit
    integer :: ierr
    integer :: n_e, n, i
    real(kind=dp) ::  qc, dos 
    real(kind=dp) ::  e_scan
    complex(kind=dp) :: e_scan_cmp
    complex(kind=dp), allocatable, dimension(:,:) :: tot, tott
    complex(kind=dp), allocatable, dimension(:,:) :: g_B, gR, gL
    complex(kind=dp), allocatable, dimension(:,:) :: sLr, sRr
    complex(kind=dp), allocatable, dimension(:,:) :: s1, s2, c1
 !   character(len=50) :: filename
    character(len=9)  :: cdate, ctime
 
    allocate (tot(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating tot in tran_bulk')
    allocate (tott(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating tott in tran_bulk')
    allocate (g_B(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating g_B in tran_bulk')
    allocate (gL(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating gL in tran_bulk')
    allocate (gR(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating gR in tran_bulk')
    allocate (sLr(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating sLr in tran_bulk')
    allocate (sRr(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating sRr in tran_bulk')
    allocate (s1(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating s1 in tran_bulk')
    allocate (s2(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating s2 in tran_bulk')
    allocate (c1(tran_num_bb,tran_num_bb),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating c1 in tran_bulk')



    !   set up the layer hamiltonians

       ! compute according to Fisher and Lee
       ! retarded Green

       e_scan_cmp = e_scan+eta 
       call tran_transfer(tot,tott,hB0,hB1,e_scan_cmp,tran_num_bb)
       call tran_green(tot,tott,hB0,hB1,e_scan,g_B,0,1,tran_num_bb) 
 
       ! compute S_Lr and S_Rr

       c1(:,:) = cmplx(hB1(:,:),kind=dp)

       ! Self-energy (Sigma_L^r) : sLr = (hB1)^+ * tott
       ! Self-energy (Sigma_R^r) : sRr = (hB1)   * tot
       sLr = cmplx_0
       sRr = cmplx_0
       call ZGEMM('C','N',tran_num_bb,tran_num_bb,tran_num_bb,cmplx_1,c1,tran_num_bb,tott,tran_num_bb,cmplx_0,sLr,tran_num_bb)
       call ZGEMM('N','N',tran_num_bb,tran_num_bb,tran_num_bb,cmplx_1,c1,tran_num_bb,tot, tran_num_bb,cmplx_0,sRr,tran_num_bb)

       ! Gamma_L = i(Sigma_L^r-Sigma_L^a)
       gL = cmplx_i*(sLr - conjg(transpose(sLr)))
       ! Gamma_R = i(Sigma_R^r-Sigma_R^a)
       gR = cmplx_i*(sRr - conjg(transpose(sRr)))

       s1 = cmplx_0
       s2 = cmplx_0
       c1 = cmplx_0
       ! s1 = Gamma_L * g_B^r
       call ZGEMM('N','N',tran_num_bb,tran_num_bb,tran_num_bb,cmplx_1,gL,tran_num_bb,g_B,tran_num_bb,cmplx_0,s1,tran_num_bb)
       ! s2 = Gamma_L * g_B^r * Gamma_R 
       call ZGEMM('N','N',tran_num_bb,tran_num_bb,tran_num_bb,cmplx_1,s1,tran_num_bb,gR, tran_num_bb,cmplx_0,s2,tran_num_bb)
       ! c1 = Gamma_L * g_B^r * Gamma_R * g_B^a
       call ZGEMM('N','C',tran_num_bb,tran_num_bb,tran_num_bb,cmplx_1,s2,tran_num_bb,g_B,tran_num_bb,cmplx_0,c1,tran_num_bb)
           
       qc = 0.0_dp
       do i=1,tran_num_bb
          qc = qc + real(c1(i,i),dp)
       end do


       dos = 0.0_dp
       do i=1,tran_num_bb
          dos = dos - aimag(g_B(i,i))
       end do
       dos = dos / pi



    deallocate (c1,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating c1 in tran_bulk')
    deallocate (s2,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating s2 in tran_bulk')
    deallocate (s1,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating s1 in tran_bulk')
    deallocate (sRr,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating sRr in tran_bulk')
    deallocate (sLr,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating sLr in tran_bulk')
    deallocate (gR,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating gR in tran_bulk')
    deallocate (gL,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating gL in tran_bulk')
    deallocate (g_B,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating g_B in tran_bulk')
    deallocate (tott,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating tott in tran_bulk')
    deallocate (tot,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating tot in tran_bulk')


    return

  end subroutine surface_dos 
  
  
  
  
   !==================================================================!
  subroutine tran_transfer(tot,tott,h_00,h_01,e_scan_cmp,nxx) 
    !==================================================================!
    !                                                                  !
    ! iterative construction of the transfer matrix                    !
    ! as Lopez-Sancho^2&Rubio, J.Phys.F:Met.Phys., v.14, 1205 (1984)   !
    ! and ibid. v.15, 851 (1985)                                       !
    !                                                                  !
    !===================================================================


    implicit none

    integer, intent(in) :: nxx
    complex(kind=dp), intent(in) ::  e_scan_cmp
    complex(kind=dp), intent(out) ::  tot(nxx,nxx)
    complex(kind=dp), intent(out) ::  tott(nxx,nxx)
    complex(kind=dp), intent(in) :: h_00(nxx,nxx)
    complex(kind=dp), intent(in) :: h_01(nxx,nxx)
    
    !
    integer  :: ierr, info
    integer  :: i, j, n, nxx2
    integer, allocatable :: ipiv(:)
    real(kind=dp) :: conver,conver2
    complex(kind=dp), allocatable, dimension(:,:) :: tsum, tsumt
    complex(kind=dp), allocatable, dimension(:,:) :: t11, t12
    complex(kind=dp), allocatable, dimension(:,:) :: s1, s2
    complex(kind=dp), allocatable, dimension(:,:,:) :: tau, taut

    allocate(ipiv(nxx),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ipiv in tran_transfer')
    allocate(tsum(nxx,nxx),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating tsum in tran_transfer')
    allocate(tsumt(nxx,nxx),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating tsumt in tran_transfer')
    allocate(t11(nxx,nxx),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating t11 in tran_transfer')
    allocate(t12(nxx,nxx),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating t12 in tran_transfer')
    allocate(s1(nxx,nxx),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating s1 in tran_transfer')
    allocate(s2(nxx,nxx),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating s2 in tran_transfer')
    allocate(tau(nxx,nxx,2),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating tau in tran_transfer')
    allocate(taut(nxx,nxx,2),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating taut in tran_transfer')

    nxx2 = nxx*nxx
    
    tot = cmplx_0
    tott = cmplx_0

    ! construction of the transfer matrix
    ! t12 = e - h_00
    t12(:,:) = cmplx(-h_00(:,:),kind=dp)
    do i=1,nxx
       t12(i,i) = e_scan_cmp + t12(i,i)
    end do

    ! compute (e - h_00)^-1 and store it in t11 
    t11 = cmplx_0
    do i=1,nxx
       t11(i,i) = cmplx_1
    end do

    ! inverse of t12 -> t11
    call ZGESV(nxx,nxx,t12,nxx,ipiv,t11,nxx,info)
    if (info.ne.0) then
       write(*,*) 'ERROR:  IN ZGESV IN tran_transfer, INFO=', info 
       call io_error('tran_transfer: problem in ZGESV 1') 
    end if

    ! compute intermediate t-matrices (defined as tau(nxx,nxx,niter)
    ! and taut(...)):
    tau = cmplx_0
    taut = cmplx_0

    ! t_0:
    t12(:,:) = cmplx(h_01(:,:),kind=dp)

    !  tau  = ( e - H_00 )^-1 * H_01^+
    call ZGEMM('N','C',nxx,nxx,nxx,cmplx_1,t11,nxx,t12,nxx,cmplx_0,tau(1,1,1),nxx)
    !  taut = ( e - H_00 )^-1 * H_01
    call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,t11,nxx,t12,nxx,cmplx_0,taut(1,1,1),nxx)

    !   initialize T:
    tot(:,:) =tau(:,:,1)
    tsum(:,:)=taut(:,:,1)

    !   initialize T^bar:
    tott(:,:) =taut(:,:,1)
    tsumt(:,:)=tau(:,:,1)

    !   main loop:
    do n=1,nterx
                
       t11 = cmplx_0
       t12 = cmplx_0

       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,tau(1,1,1),nxx,taut(1,1,1),nxx,cmplx_0,t11,nxx)
       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,taut(1,1,1),nxx,tau(1,1,1),nxx,cmplx_0,t12,nxx)

       s1(:,:) = -t11(:,:)-t12(:,:)
       do i=1,nxx
          s1(i,i) = cmplx_1+s1(i,i)
       end do

       s2 = cmplx_0
       do i=1,nxx
          s2(i,i)=cmplx_1 
       end do

       call ZGESV(nxx,nxx,s1,nxx,ipiv,s2,nxx,info)
       if (info.ne.0) then
          write(*,*) 'ERROR:  IN ZGESV IN tran_transfer, INFO=', info 
          call io_error('tran_transfer: problem in ZGESV 2') 
       end if

       t11 = cmplx_0
       t12 = cmplx_0

       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,tau(1,1,1),nxx,tau(1,1,1),nxx,cmplx_0,t11,nxx)
       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,taut(1,1,1),nxx,taut(1,1,1),nxx,cmplx_0,t12,nxx)
       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,s2,nxx,t11,nxx,cmplx_0,tau(1,1,2),nxx)
       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,s2,nxx,t12,nxx,cmplx_0,taut(1,1,2),nxx)

       !   put the transfer matrices together

       t11 = cmplx_0
       s1  = cmplx_0

       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,tsum,nxx,tau(1,1,2),nxx,cmplx_0,t11,nxx)
       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,tsum,nxx,taut(1,1,2),nxx,cmplx_0,s1,nxx)
       call ZCOPY(nxx2,t11,1,s2,1)
       call ZAXPY(nxx2,cmplx_1,tot,1,s2,1)

       tot(:,:) = s2(:,:)
       tsum(:,:)= s1(:,:)

       t11 = cmplx_0
       s1  = cmplx_0

       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,tsumt,nxx,taut(1,1,2),nxx,cmplx_0,t11,nxx)
       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,tsumt,nxx,tau(1,1,2),nxx,cmplx_0,s1,nxx) 
       call ZCOPY(nxx2,t11,1,s2,1)
       call ZAXPY(nxx2,cmplx_1,tott,1,s2,1)

       tott(:,:) = s2(:,:)
       tsumt(:,:)= s1(:,:)

       tau(:,:,1) = tau(:,:,2)
       taut(:,:,1)= taut(:,:,2)

       ! convergence check on the t-matrices

       conver = 0.0_dp
       conver2 = 0.0_dp

       do j=1,nxx
          do i=1,nxx
              conver=conver+sqrt(real(tau(i,j,2),dp)**2+aimag(tau(i,j,2))**2)
              conver2=conver2+sqrt(real(taut(i,j,2),dp)**2+aimag(taut(i,j,2))**2)
          end do
       end do

       if (conver.lt.eps7 .and. conver2.lt.eps7) return
    end do 

    if (conver.gt.eps7 .or. conver2.gt.eps7) &
       call io_error('Error in converging transfer matrix in tran_transfer') 

    deallocate(ipiv,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating ipiv in tran_transfer')
    deallocate(tsum,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating tsum in tran_transfer')
    deallocate(tsumt,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating tsumt in tran_transfer')
    deallocate(t11,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating t11 in tran_transfer')
    deallocate(t12,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating t12 in tran_transfer')
    deallocate(s1,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating s1 in tran_transfer')
    deallocate(s2,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating s2 in tran_transfer')
    deallocate(tau,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating tau in tran_transfer')
    deallocate(taut,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating taut in tran_transfer')

    return 

  end subroutine tran_transfer

  !==================================================================!
  subroutine tran_green(tot,tott,h_00,h_01,e_scan, g,igreen,invert,nxx)
    !==================================================================!
    !   construct green's functions
    !   
    !   igreen = -1  left surface
    !   igreen =  1  right surface
    !   igreen =  0  bulk
   
    !   invert = 0 computes g^-1
    !   invert = 1 computes g^-1 and g
    !==================================================================!



    implicit none

    integer, intent(in) :: nxx
    integer, intent(in) :: igreen
    integer, intent(in) :: invert
    real(kind=dp),     intent(in) :: e_scan 

    complex(kind=dp),  intent(in) :: tot(nxx,nxx),tott(nxx,nxx)
    complex(kind=dp), intent(out) :: g(nxx,nxx)

    complex(kind=dp),     intent(in) :: h_00(nxx,nxx), h_01(nxx,nxx)


    integer :: ierr, info
    integer :: i
    integer, allocatable :: ipiv(:)
    complex(kind=dp), allocatable, dimension(:,:) :: g_inv, eh_00
    complex(kind=dp), allocatable, dimension(:,:) :: s1, s2, c1           

    allocate(ipiv(nxx),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ipiv in tran_green')
    allocate(g_inv(nxx,nxx))
    if (ierr/=0) call io_error('Error in allocating g_inv in tran_green')
    allocate(eh_00(nxx,nxx))
    if (ierr/=0) call io_error('Error in allocating eh_00 in tran_green')
    allocate(c1(nxx,nxx))
    if (ierr/=0) call io_error('Error in allocating c1 in tran_green')
    allocate(s1(nxx,nxx))
    if (ierr/=0) call io_error('Error in allocating s1 in tran_green')
    allocate(s2(nxx,nxx))
    if (ierr/=0) call io_error('Error in allocating s2 in tran_green')

    c1(:,:)=cmplx(h_01(:,:),kind=dp)

    select case(igreen)

      case(1) 

       ! construct the surface green's function g00 

       s1 = cmplx_0
       ! s1 = H_01 * T
       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,c1,nxx,tot,nxx,cmplx_0,s1,nxx)

       ! eh_00 =  -H_00 - H_01*T
       eh_00(:,:) = cmplx(-h_00(:,:),kind=dp)-s1(:,:)
       ! eh_00 = e_scan -H_00 - H_01*T
       do i=1,nxx
          eh_00(i,i) = cmplx(e_scan,kind=dp) + eh_00(i,i)
       end do

       g_inv(:,:) = eh_00(:,:)

       ! identity
       g = cmplx_0
       do i=1,nxx
          g(i,i)=cmplx_1 
       end do
    
       if (invert.eq.1) then
          call ZGESV(nxx,nxx,eh_00,nxx,ipiv,g,nxx,info)
          if (info.ne.0) then
             write(*,*) 'ERROR:  IN ZGESV IN tran_green, INFO=', info 
             call io_error('tran_green: problem in ZGESV 1') 
          end if
       end if

      case(-1)

       !  construct the dual surface green's function gbar00 

       s1 = cmplx_0
       ! s1 = H_01^+ * T^bar
       call ZGEMM('C','N',nxx,nxx,nxx,cmplx_1,c1,nxx,tott,nxx,cmplx_0,s1,nxx)

       ! s1 = -H_00 - H_01^+ * T^bar
       eh_00(:,:) = cmplx(-h_00(:,:),kind=dp)-s1(:,:)
       ! s1 = e_scan - H_00 - H_01^+ * T^bar
       do i=1,nxx
          eh_00(i,i) = cmplx(e_scan,kind=dp) + eh_00(i,i)
       end do
     
       g_inv(:,:) = eh_00(:,:)

       ! identity
       g = cmplx_0
       do i=1,nxx
          g(i,i)=cmplx_1
       end do

       if (invert.eq.1) then
          call ZGESV(nxx,nxx,eh_00,nxx,ipiv,g,nxx,info)
          if (info.ne.0) then
             write(*,*) 'ERROR:  IN ZGESV IN tran_green, INFO=', info 
             call io_error('tran_green: problem in ZGESV 2') 
          end if
       end if
    
      case(0)

      !  construct the bulk green's function gnn or (if surface=.true.) the
      !  sub-surface green's function

       s1 = cmplx_0
       s2 = cmplx_0
       ! s1 = H_01 * T
       call ZGEMM('N','N',nxx,nxx,nxx,cmplx_1,c1,nxx,tot,nxx,cmplx_0,s1,nxx)
       ! s2 = H_01^+ * T^bar
       call ZGEMM('C','N',nxx,nxx,nxx,cmplx_1,c1,nxx,tott,nxx,cmplx_0,s2,nxx)

       eh_00(:,:) = cmplx(-h_00(:,:),kind=dp)-s1(:,:)-s2(:,:)
       do i=1,nxx
          eh_00(i,i) = cmplx(e_scan,kind=dp) + eh_00(i,i)
       end do

       g_inv(:,:) = eh_00(:,:)

       ! identity
       g = cmplx_0
       do i=1,nxx
          g(i,i)=cmplx_1 
       end do
     
       if (invert.eq.1) then
          call ZGESV(nxx,nxx,eh_00,nxx,ipiv,g,nxx,info)
          if (info.ne.0) then
             write(*,*) 'ERROR:  IN ZGESV IN tran_green, INFO=', info 
             call io_error('tran_green: problem in ZGESV 3') 
          end if
       end if

    end select

    deallocate(s2)
    if (ierr/=0) call io_error('Error in deallocating s2 in tran_green')
    deallocate(s1)
    if (ierr/=0) call io_error('Error in deallocating s1 in tran_green')
    deallocate(c1)
    if (ierr/=0) call io_error('Error in deallocating c1 in tran_green')
    deallocate(eh_00)
    if (ierr/=0) call io_error('Error in deallocating eh_00 in tran_green')
    deallocate(g_inv)
    if (ierr/=0) call io_error('Error in deallocating g_inv in tran_green')
    deallocate(ipiv)
    if (ierr/=0) call io_error('Error in deallocating ipiv in tran_green')

    return

  end subroutine tran_green  
 
 

 

 
 
 
    !===================================================================
         subroutine utility_recip_lattice (real_lat,recip_lat,volume)  !
    !==================================================================!
    !                                                                  !
    !  Calculates the reciprical lattice vectors and the cell volume   !
    !                                                                  !
    !===================================================================


    implicit none
    real(kind=dp), intent(in)  :: real_lat (3, 3)
    real(kind=dp), intent(out) :: recip_lat (3, 3)  
    real(kind=dp), intent(out) :: volume

    recip_lat(1,1)=real_lat(2,2)*real_lat(3,3)-real_lat(3,2)*real_lat(2,3)
    recip_lat(1,2)=real_lat(2,3)*real_lat(3,1)-real_lat(3,3)*real_lat(2,1)
    recip_lat(1,3)=real_lat(2,1)*real_lat(3,2)-real_lat(3,1)*real_lat(2,2)
    recip_lat(2,1)=real_lat(3,2)*real_lat(1,3)-real_lat(1,2)*real_lat(3,3)
    recip_lat(2,2)=real_lat(3,3)*real_lat(1,1)-real_lat(1,3)*real_lat(3,1)
    recip_lat(2,3)=real_lat(3,1)*real_lat(1,2)-real_lat(1,1)*real_lat(3,2)
    recip_lat(3,1)=real_lat(1,2)*real_lat(2,3)-real_lat(2,2)*real_lat(1,3)
    recip_lat(3,2)=real_lat(1,3)*real_lat(2,1)-real_lat(2,3)*real_lat(1,1)
    recip_lat(3,3)=real_lat(1,1)*real_lat(2,2)-real_lat(2,1)*real_lat(1,2)

    volume=real_lat(1,1)*recip_lat(1,1) + &
         real_lat(1,2)*recip_lat(1,2) + &
         real_lat(1,3)*recip_lat(1,3)  


    if( abs(volume) < eps5 ) then
       call io_error(' Found almost zero Volume in utility_recip_lattice')
    end if

    recip_lat=twopi*recip_lat/volume
    volume=abs(volume)

    return

         end subroutine utility_recip_lattice
         
         
         
         
         

    !===================================================================
         subroutine utility_compute_metric(real_lat,recip_lat, &
                                              real_metric,recip_metric)
    !==================================================================!
    !                                                                  !
    !  Calculate the real and reciprical space metrics                 !
    !                                                                  !
    !===================================================================  
    implicit none

    real(kind=dp), intent(in)  :: real_lat(3,3)
    real(kind=dp), intent(in)  :: recip_lat(3,3)
    real(kind=dp), intent(out) :: real_metric(3,3)
    real(kind=dp), intent(out) :: recip_metric(3,3)

    integer :: i,j,l

    real_metric=0.0_dp ; recip_metric=0.0_dp

    do j=1,3
       do i=1,j
          do l=1,3
             real_metric(i,j)=real_metric(i,j)+real_lat(i,l)*real_lat(j,l)
             recip_metric(i,j)=recip_metric(i,j)+recip_lat(i,l)*recip_lat(j,l)
          enddo
          if(i.lt.j) then
             real_metric(j,i)=real_metric(i,j)
             recip_metric(j,i)=recip_metric(i,j)
          endif
       enddo
    enddo

         end subroutine utility_compute_metric
         
    !===================================================================
         subroutine utility_frac_to_cart(frac,cart,real_lat)
    !==================================================================!
    !                                                                  !
    !  Convert from fractional to Cartesian coordinates                !
    !                                                                  !
    !===================================================================  
    implicit none

    real(kind=dp), intent(in)  :: real_lat(3,3)
    real(kind=dp), intent(in)  :: frac(3)
    real(kind=dp), intent(out) :: cart(3)

    integer :: i

    do i=1,3
       cart(i)=real_lat(1,i)*frac(1) + real_lat(2,i)*frac(2) + real_lat(3,i)*frac(3) 
    end do

    return

  end subroutine utility_frac_to_cart


    !===================================================================
         subroutine utility_cart_to_frac(cart,frac,recip_lat)
    !==================================================================!
    !                                                                  !
    !  Convert from fractional to Cartesian coordinates                !
    !                                                                  !
    !===================================================================  

    implicit none

    real(kind=dp), intent(in)  :: recip_lat(3,3)
    real(kind=dp), intent(out)  :: frac(3)
    real(kind=dp), intent(in)  :: cart(3)

    integer :: i

    do i=1,3
       frac(i)=recip_lat(i,1)*cart(1) + recip_lat(i,2)*cart(2) + recip_lat(i,3)*cart(3) 
    end do

    frac=frac/twopi


    return

  end subroutine utility_cart_to_frac         
         
         
    !========================================================!
    subroutine utility_translate_home(vec,real_lat,recip_lat)
    !========================================================!
    !                                                        !
    !        Translate a vector to the home unit cell        !
    !                                                        !
    !========================================================!

      implicit none

      real(kind=dp), intent(inout) :: vec(3)
      real(kind=dp), intent(in)    :: real_lat(3,3)
      real(kind=dp), intent(in)    :: recip_lat(3,3)
      
      ! <<<local variables>>>
      integer       :: ind
      real(kind=dp) :: r_home(3),r_frac(3)
      real(kind=dp) :: shift

      r_home=0.0_dp;r_frac=0.0_dp

      ! Cartesian --> fractional
      call utility_cart_to_frac(vec,r_frac,recip_lat)
      ! Rationalise to interval [0,1]
      do ind=1,3
         if (r_frac(ind).lt.0.0_dp) then
            shift=real(ceiling(abs(r_frac(ind))),kind=dp)
            r_frac(ind)=r_frac(ind)+shift
         endif
         if (r_frac(ind).gt.1.0_dp) then
            shift=-real(int(r_frac(ind)),kind=dp)
            r_frac(ind)=r_frac(ind)+shift
         endif
      enddo
      ! Fractional --> Cartesian
      call utility_frac_to_cart(r_frac,r_home,real_lat)
      
      vec = r_home

      return
    end subroutine utility_translate_home       


    end program surface_stat
