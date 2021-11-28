program lorentzian_for_omet

    implicit none

    CHARACTER (len=256) :: outdir                                                       ! Current directry of output
    CHARACTER (len=256) :: outdir_matele_opt, outdir_eigenvalues                        ! Directory for matrix element transition and eigenvalues
    CHARACTER (len=256) :: filename_matele_opt                                          ! Optical matrix element transition file use for calculation
    CHARACTER (len=256) :: valence_bands, conduction_bands                              ! Valence and conduction eigenvalues file
    CHARACTER (len=256) :: n1, n2                                                       ! Band number to character for reading the eigenvalues file

    integer :: ik, n, z, setk
    integer :: nec, nev                                                                 ! Valence and conduction band number

    real (KIND=KIND(1.0d0)) :: el1, el2, el3, ef                                        ! Input for 3 laser energy and fermi energy
    real (KIND=KIND(1.0d0)) :: gamma                                                    ! Level broadening in eV
    real (KIND=KIND(1.0d0)) :: pi                                                       ! Pi number
    real (KIND=KIND(1.0d0)) :: ng, ngp                                                  ! Simplified constant

    real , dimension(1000) :: kx, ky, kz, abslft, absrh                                 ! Optical matrix element transition parameter
    real , dimension(1000) :: x1, y1, z1, ev                                            ! Valence band parameter
    real , dimension(1000) :: x2, y2, z2, ec                                            ! Conduction band parameter 

    ! Initial message
    write(0,*) ''
    write(0,*) '---------------------------------------------------------------------'
    write(0,*) 'Lorentzian function as delta function'
    write(0,*) 'Created by R. Layton & MJ. Prakoso, started from 2020/09/04'
    write(0,*) '---------------------------------------------------------------------'
    write(0,*) ''
    write(0,*) 'There are 3 laser energy in eV'
    write(0,*) 'Please input the needed parameter that needed for calculation'
    write(0,*) '* symbol represent the needed parameter for calculation'
    write(0,*) ''
    write(0,*) '---------------------------------------------------------------------'
    write(0,*) ''

    ! Read laser energy
    write(0,*) 'Please input the 1st laser energy in eV'
    read (*,'(e16.8e3)') el1
    write(0,*) ''
    write(21,*) '# 1st Laser energy selected : ', el1

    write(0,*) 'Please input the 2nd laser energy in eV'
    read (*,'(e16.8e3)') el2
    write(0,*) ''
    write(21,*) '# 2nd Laser energy selected : ', el2

    write(0,*) 'Please input the 3rd laser energy in eV'
    read (*,'(e16.8e3)') el3
    write(0,*) ''
    write(21,*) '# 3rd Laser energy selected : ', el3
    write(21,*) ''

    ! Read fermi energy
    write(0,*) '*Please input the Fermi level energy in eV'
    read (*,'(e16.8e3)') ef
    write(0,*) ''
    write(21,*) '# Fermi level energy selected : ', ef
    write(21,*) ''

    ! Read gamma parameter
    write(0,*) '*Please input the level broadening in eV for Lorentzian Function'
    read (*,'(e16.8e3)') gamma
    write(0,*) ''
    write(21,*) '# Level broadening selected : ', gamma

    ! Read valence band number
    write(0,*) '*Please input the valence band number'
    read *, nev
    write(0,*) ''

    ! Read conduction band number
    write(0,*) '*Please input the conduction band number'
    read *, nec
    write(0,*) ''

    write(21,*) '# Valence band number : ', nev, ', Conduction band number : ', nec
    write(21,*) ''

    ! Read grid for intensity 
    write(0,*) ''
    write(0,*) 'For plotting intensity of polarization'
    write(0,*) 'This features is ONLY for square grid with hexagonal lattice'
    write(0,*) ''
    write(0,*) 'Example: If grid is 30x30x1 input the k-point grid as 30'
    write(0,*) 'Please input the grid for calculating the intensity'
    write(0,*) 'Input 0 for the grid input to abandon this features'
    read *, setk
    write(0,*) ''
    write(21,*) '# K-point grid : ', setk

    ! Reading matrix element transition and eigenvalues directories and files 
    outdir = '.'
    WRITE(n1,'(I3)') nev
    WRITE(n2,'(I3)') nec

    outdir_eigenvalues = TRIM(ADJUSTL(outdir))//'/eigenvalues/'
    valence_bands = TRIM(ADJUSTL(outdir_eigenvalues))//'eigenvalues_'//TRIM(ADJUSTL(n1))//'.dat'
    conduction_bands = TRIM(ADJUSTL(outdir_eigenvalues))//'eigenvalues_'//TRIM(ADJUSTL(n2))//'.dat'

    outdir_matele_opt = TRIM(ADJUSTL(outdir))//'/matele_opt/'
    filename_matele_opt = TRIM(ADJUSTL(outdir_matele_opt))//'matele_opt_'//TRIM(ADJUSTL(n1))//'_'//TRIM(ADJUSTL(n2))//'.dat'

    OPEN(12, file=filename_matele_opt)
    OPEN(15, file=valence_bands)
    OPEN(18, file=conduction_bands)

    ! Calculation for lorentzian function
    OPEN(21, file='opt_met.dat')

    ! Definition of Lorentzian = ((0.5*gamma)/pi) / (((EC - EV)**2) + ((0.5*gamma)**2))
    ! EC and EV as 1st conduction and valence band eigenvalues
    ! Simplify with (0.5*gamma) as ng, ((0.5*gamma)/pi) as ngp

    write(0,*) 'Starting calculation'
    write(0,*) ''

    pi = 3.14156
    ng = (0.5*gamma)
    ngp = ng/pi
    write(21,*) '# kx ky kz Left(elaser1) Right(elaser1) Left(elaser2) Right(elaser2) Left(elaser3) Right(elaser3)'

    do ik = 1,10000
        read(12,'(3(1x,f10.4),6(1x,e16.8e3))') kx(ik), ky(ik), kz(ik), abslft(ik), absrh(ik)
        read(15,'(3(1x,f10.4),1(1x,e16.8e3))') x1(ik), y1(ik), z1(ik), ev(ik)
        read(18,'(3(1x,f10.4),1(1x,e16.8e3))') x2(ik), y2(ik), z2(ik), ec(ik)

        if (ev(ik) >= ef) then
            ev(ik) = ef
        else if (ef == 0) then
            ev(ik) = ev(ik)
        else if (ev(ik) <= ef) then
            ev(ik) = ev(ik)
        end if

        write(21,'(3(1x,f10.4),6(1x,e16.8e3))') kx(ik), ky(ik), kz(ik), &
                abslft(ik)*ngp/((((ec(ik)-ev(ik))-el1)**2)+(ng**2)), &
                absrh(ik)*ngp/((((ec(ik)-ev(ik))-el1)**2)+(ng**2)), &
                abslft(ik)*ngp/((((ec(ik)-ev(ik))-el2)**2)+(ng**2)), &
                absrh(ik)*ngp/((((ec(ik)-ev(ik))-el2)**2)+(ng**2)), &
                abslft(ik)*ngp/((((ec(ik)-ev(ik))-el3)**2)+(ng**2)), &
                absrh(ik)*ngp/((((ec(ik)-ev(ik))-el3)**2)+(ng**2))

        if (setk == 0) then
            goto 4
        else
            goto 3
        end if

        ! Extract data for intensity plot
3       OPEN(25, file='degree_pol.dat')
        if (ik == 1) then
            write(25,*) '# kx, ky, Ev, Ec, Left(elaser1) Right(elaser1) Left(elaser2) Right(elaser2) Left(elaser3) Right(elaser3)'
        end if

        do z = 0,setk
            n = 1 + setk*z + z
            if (ik == n) then
                write(25,'(4(1x,f10.4),9(1x,e16.8e3))') kx(ik), ky(ik), ev(ik), ec(ik), &
                    abslft(ik)*ngp/((((ec(ik)-ev(ik))-el1)**2)+(ng**2)), &
                    absrh(ik)*ngp/((((ec(ik)-ev(ik))-el1)**2)+(ng**2)), &
                    abslft(ik)*ngp/((((ec(ik)-ev(ik))-el2)**2)+(ng**2)), &
                    absrh(ik)*ngp/((((ec(ik)-ev(ik))-el2)**2)+(ng**2)), &
                    abslft(ik)*ngp/((((ec(ik)-ev(ik))-el3)**2)+(ng**2)), &
                    absrh(ik)*ngp/((((ec(ik)-ev(ik))-el3)**2)+(ng**2))
            end if
        end do
4   end do

    write(0,*) 'Calculation done'

end program lorentzian_for_omet