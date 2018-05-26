!************************************************************
! This program calculated 3d Fourier transform for arbitray *
! dimensions using FFTW library BASIC interface-            *
! Fortran interface & C interface.                          *
! While using C interface dimension are reversed accounting *
! to difference of C & Fortran arrays.                      *
!************************************************************

        module fftw_basic_interface
        contains

        subroutine fftw_3d(input)

#ifdef __CINTERFACE__
        use, intrinsic :: iso_c_binding
#endif

        implicit none
        integer nx,ny,nz
        parameter(nx=4,ny=4,nz=4)

#ifdef __FINTERFACE__
        integer*8 plan_3d_forward, plan_3d_backward
        double complex input(nx,ny,nz)
#endif

#ifdef __CINTERFACE__
        type(C_PTR)  plan_3d_forward, plan_3d_backward
        complex(C_DOUBLE_COMPLEX)  input(nx,ny,nz)
#endif

#ifdef __CINTERFACE__
        include 'fftw3.f03'
#endif

#ifdef __FINTERFACE__
        include 'fftw3.f'
#endif

! Forward 3 D FFT
#ifdef __FINTERFACE__
        call dfftw_plan_dft_3d(plan_3d_forward,nx,ny,nz,input,input,&
                  FFTW_FORWARD,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan_3d_forward,input,input)
        call dfftw_destroy_plan(plan_3d_forward)
#endif

#ifdef __CINTERFACE__
        plan_3D_forward = fftw_plan_dft_3d(nz,ny,nx,input,input, &
                  FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(plan_3d_forward,input,input)
        call fftw_destroy_plan(plan_3d_forward)
#endif



! Backward 3D FFT
#ifdef __FINTERFACE__
        call dfftw_plan_dft_3d(plan_3d_backward,nx,ny,nz,input,input,&
                  FFTW_BACKWARD,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan_3d_backward,input,input)
        call dfftw_destroy_plan(plan_3d_backward)
#endif

#ifdef __CINTERFACE__
        plan_3d_backward = fftw_plan_dft_3d(nz,ny,nx,input,input,&
                  FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(plan_3d_backward,input,input)
        call fftw_destroy_plan(plan_3d_backward)
#endif

        end
        end 

        program ft
        use fftw_basic_interface
        implicit none 
        include 'fftw3.f'

        integer n,i
        parameter(N=64)
        double complex input(N)

! Data filling 

        do i=1,n
        input(i) = i  
        end do

        print*, 'Data before forward fft using fftw fortran interface',&
                input


        print*, 'Data before forward fft using fftw C interface',input


        call fftw_3d(input)

! Normalization -Note -FFTW produces un-normalized output
        do i=1,n
           input(i) = input(i)/DBLE(N)
        end do

#ifdef __FINTERFACE__
        print*, 'Data after backward fft using fftw BASIC fortran'&
                 'interface', input
#endif

#ifdef __CINTERFACE__
        print*, 'Data after backward fft using fftw BASIC C interface'&
                 , input
#endif

        end
