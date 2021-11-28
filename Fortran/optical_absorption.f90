!----------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE inputoutput_mod
!----------------------------------------------------------------------------------------------------------------------------------------------------------
    
    IMPLICIT NONE
        
    PRIVATE
    SAVE
        
    PUBLIC :: stdout, &
        xml_pw_export, &
        fn_matele_opt, fn_eigv
        
    ! Output of standard output
    INTEGER :: stdout = 6                                               ! Standard output

    ! Reading data file
    INTEGER :: xml_pw_export = 102                                      ! Xml file created by pw_export in Quantum Espresso

    ! Output of data files
    INTEGER :: fn_matele_opt = 203                                      ! Output for the electron-photon matrix element
    INTEGER :: fn_eigv = 206                                            ! Output for eigen energy list
    
!----------------------------------------------------------------------------------------------------------------------------------------------------------
END MODULE inputoutput_mod
!----------------------------------------------------------------------------------------------------------------------------------------------------------
    
!----------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE dipolevector_mod
!----------------------------------------------------------------------------------------------------------------------------------------------------------
        
    IMPLICIT NONE
        
    PRIVATE
    SAVE
        
    PUBLIC :: unit_parameter, polarization_vector, &
            prefix, outdir, input_pw_export, &
            sorb, circular_pol, nonpol, plot_matele_opt, &
            valence_band, conduction_band, norm, &
            nks, nbnd, ngs, nspin, npol, &
            alat, &
            k, g, wk, eigv, avec, bvec, &
            dvec, matele_opt, coef_wfc, polvec, &
            pi, im, &
            ry2ev, ev2ry, ev2cm, ry2cm
        
    CHARACTER (len=256) :: prefix                                       ! Prefix of the input/output data
    CHARACTER (len=256) :: outdir                                       ! Directry of output
    CHARACTER (len=256) :: input_pw_export                              ! File name of pw_export xml file
    
    LOGICAL :: sorb                                                     ! Consideration of spin-ortbit interaction (sorb = .FALSE. in this version)
    LOGICAL :: circular_pol                                             ! If true, circular polarizaed incident light is calculated
    LOGICAL :: nonpol                                                   ! If true, calculate for non-polarized light
    LOGICAL :: plot_matele_opt                                          ! If true, generate the file of matrix element for electron-photon matrix element

    INTEGER :: valence_band                                             ! Valence band number for dipole vector calculation
    INTEGER :: conduction_band                                          ! Conduction band number for dipole vector calculation
    INTEGER :: norm                                                     ! Normalization length and size factor for dipole vector calculation
    
    INTEGER :: nks                                                      ! Total number of k points in scf calculation by Quantum Espresso
    INTEGER :: nbnd                                                     ! Total number of energy bands
    INTEGER :: ngs                                                      ! Total number of plane wave for wave functions
    INTEGER :: nspin                                                    ! Number of spin index (1:without 2:with spin orbit interaction)
    INTEGER :: npol                                                     ! Number of polarization vector
        
    REAL (KIND=KIND(1.0d0)) :: alat                                     ! Lattice constant (atomic unit, Bohr)
    
    REAL (KIND=KIND(1.0d0)), ALLOCATABLE :: k(:,:)                      ! Wave vector k(ik,idim)
    REAL (KIND=KIND(1.0d0)), ALLOCATABLE :: g(:,:)                      ! Reciprocal lattice vector g(ig,idim)
    REAL (KIND=KIND(1.0d0)), ALLOCATABLE :: wk(:)                       ! Weight of each k point wk(ik)
    REAL (KIND=KIND(1.0d0)), ALLOCATABLE :: eigv(:,:)                   ! Eigen energy of each band at each k point eigv(ik,ib)
    REAL (KIND=KIND(1.0d0)), ALLOCATABLE :: avec(:,:)                   ! Primitive lattice vector
    REAL (KIND=KIND(1.0d0)), ALLOCATABLE :: bvec(:,:)                   ! Primitive reciplocal lattice vector
        
    COMPLEX (KIND=KIND(1.0d0)), ALLOCATABLE :: dvec(:,:,:,:)            ! Transition dipole moment (dipole vector) dvec(ik,ib_i,ib_f,ixyz)
    COMPLEX (KIND=KIND(1.0d0)), ALLOCATABLE :: matele_opt(:,:,:,:)      ! Optical matrix element matele_opt(ik,ibi,ibf,ip)
    COMPLEX (KIND=KIND(1.0d0)), ALLOCATABLE :: coef_wfc(:,:,:,:)        ! Coefficient of wave function coef_wfc(ik,ib,ig,ispin)
    COMPLEX (KIND=KIND(1.0d0)), ALLOCATABLE :: polvec(:,:)              ! Polarization (Jones') vector polvec(npol,3) 1:x, 2:y, 3:z, 4:LCP, 5:RCP
    
    ! parameter
    REAL (KIND=KIND(1.0d0)) pi                                          ! Circumference ratio
    COMPLEX (KIND=KIND(1.0d0)) im                                       ! Imaginary number i
    
    ! conversion of unit
    REAL (KIND=KIND(1.0d0)) ry2ev                                       ! Ry --> eV
    REAL (KIND=KIND(1.0d0)) ev2ry                                       ! eV --> Ry
    REAL (KIND=KIND(1.0d0)) ev2cm                                       ! eV --> cm-1
    REAL (KIND=KIND(1.0d0)) ry2cm                                       ! Ry --> cm-1
        
    CONTAINS
    
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE unit_parameter
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    
        pi = acos(-1.0d0)
        ry2eV = 13.60568d0
        ev2ry = 1.0d0 / ry2ev
        ev2cm = 8065.0d0
        ry2cm = ry2ev * ev2cm
            
        im = CMPLX(0.0d0,1.0d0)
    
        !WRITE(0,*) 'Conversion Done'
        
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE unit_parameter
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE polarization_vector
    !------------------------------------------------------------------------------------------------------------------------------------------------------
            
        IF (circular_pol .EQ. .FALSE. .AND.  nonpol .EQ. .FALSE.) THEN
            npol = 3

            !WRITE(0,*) 'Linear Polarized Selected, npol = 3'
        ELSEIF (circular_pol .EQ. .TRUE. .AND. nonpol .EQ. .FALSE.) THEN
            npol = 5
    
            !WRITE(0,*) 'Circular Polarized Selected, npol = 5'
        ELSEIF (nonpol .EQ. .TRUE.) THEN
            npol = 6
    
            !WRITE(0,*) 'Non Polarized Selected, npol = 6'
        ENDIF
            
        ALLOCATE(polvec(npol,3))
        polvec(:,:) = 0.0d0

        !ips = 1, p=(1,0,0)
        polvec(1,1) = cmplx(1.0d0, 0.0d0)
        polvec(1,2) = cmplx(0.0d0, 0.0d0)
        polvec(1,3) = cmplx(0.0d0, 0.0d0)
    
        !ips = 2, p=(0,1,0)
        polvec(2,1) = cmplx(0.0d0, 0.0d0)                               ! polvec(2,1) = cmplx(cos(pi/3.0d0), 0.0d0)
        polvec(2,2) = cmplx(1.0d0, 0.0d0)                               ! polvec(2,2) = cmplx(sin(pi/3.0d0), 0.0d0)
        polvec(2,3) = cmplx(0.0d0, 0.0d0)
    
        !ips = 3, p=(0,0,1)
        polvec(3,1) = cmplx(0.0d0, 0.0d0)
        polvec(3,2) = cmplx(0.0d0, 0.0d0)
        polvec(3,3) = CMPLX(1.0d0, 0.0d0)
    
        !WRITE(0,*) 'Linear Polarized Light Clear'
    
        IF (circular_pol .eq. .TRUE.) THEN
            ! ips = 4, p=(1,i,0) Left handed circular polarizaed light
            polvec(4,1) = (1.0d0/SQRT(2.0d0))*CMPLX(1.0d0, 0.0d0)
            polvec(4,2) = (1.0d0/SQRT(2.0d0))*CMPLX(0.0d0, 1.0d0)
            polvec(4,3) = CMPLX(0.0d0, 0.0d0)
            
            ! ips = 5, p=(1,-i,0) Right-handed circular polarized light
            polvec(5,1) = (1.0d0/SQRT(2.0d0))*CMPLX(1.0d0, 0.0d0)
            polvec(5,2) = (1.0d0/SQRT(2.0d0))*CMPLX(0.0d0, -1.0d0)
            polvec(5,3) = CMPLX(0.0d0, 0.0d0)
        ENDIF
    
        !WRITE(0,*) 'Circular Polarized Light Clear'
    
        IF (nonpol .EQ. .TRUE.) THEN
            polvec(6,1) = (1.0d0/SQRT(3.0d0))*CMPLX(1.0d0, 0.0d0)
            polvec(6,2) = (1.0d0/SQRT(3.0d0))*CMPLX(1.0d0, 0.0d0)
            polvec(6,3) = (1.0d0/SQRT(3.0d0))*CMPLX(1.0d0, 0.0d0)
        ENDIF
    
        !WRITE(0,*) 'Non Polarized Light Clear'
        !WRITE(0,*) 'Matrix Polarized Vector Program Clear'
            
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE polarization_vector
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    
!----------------------------------------------------------------------------------------------------------------------------------------------------------
END MODULE dipolevector_mod
!----------------------------------------------------------------------------------------------------------------------------------------------------------
    
!----------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE readingfile_mod
!----------------------------------------------------------------------------------------------------------------------------------------------------------
    
    USE  inputoutput_mod, ONLY : stdout, xml_pw_export
    USE  dipolevector_mod, ONLY : prefix, outdir, input_pw_export, &
                                sorb, circular_pol, nonpol, plot_matele_opt, &
                                nks, nbnd, ngs, nspin, &
                                valence_band, conduction_band, norm, alat, &
                                k, g, wk, eigv, avec, bvec, &
                                dvec, coef_wfc, &
                                pi, im, &
                                ev2ry
    
    IMPLICIT NONE
    
    CONTAINS
    
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE reading_QEinputfile
    !------------------------------------------------------------------------------------------------------------------------------------------------------
       
    IMPLICIT NONE
        
    INTEGER :: ios
    
    namelist / inputdvec / &
        outdir, prefix, input_pw_export, &
        sorb, circular_pol, nonpol, plot_matele_opt, &
        valence_band, conduction_band, norm
    
    ! Set default values for variables in namelist
    outdir = '.'
    prefix = ''
    input_pw_export = ''
    sorb = .FALSE.
    circular_pol = .TRUE.
    nonpol = .FALSE.
    plot_matele_opt = .TRUE.
    valence_band = 1
    conduction_band = 2
    norm = 1
    
    ! Read the namelist inputdvec
    READ (5, inputdvec, iostat = ios)
    
    !WRITE(0,*) 'Reading Input Program Done'
        
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE reading_QEinputfile
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE reading_xmlfile
    !------------------------------------------------------------------------------------------------------------------------------------------------------
            
        INTEGER :: i, u, v, w, n, nn, ispin
        INTEGER, ALLOCATABLE :: gik(:,:), n_gik(:)

        REAL(KIND=KIND(1.0d0)) :: r, s, sum_wk
        REAL(KIND=KIND(1.0d0)), ALLOCATABLE :: gi(:,:), g2(:,:)
        CHARACTER(LEN=256) :: bb, cc, dd, ee
            
        WRITE(stdout,'(1x)')
        WRITE(stdout,'(5x,"Reading xml data file")')
        WRITE(stdout,'(1x)')
            
        ! Open the pw_export.xml file
        OPEN(xml_pw_export, file=TRIM(ADJUSTL(input_pw_export)))
    
        IF (sorb == .TRUE.) then
            nspin = 2
    
            !WRITE(0,*) 'Spin Orbit Coupling Selected'
        ELSEIF (sorb == .FALSE.) then
            nspin = 1
    
            !WRITE(0,*) 'Non Spin Orbit Coupling Selected'
        ENDIF
    
        ispin = 1
    
        !WRITE(0,*) 'Number Spin =', ispin
            
        ! read nks
        n=0
        DO u = 1, 100
            READ(xml_pw_export,'(a)') bb
            IF (INDEX(bb, 'nktot="') /= 0) THEN
            w = INDEX(bb, 'nktot="') + 7
                v = w + INDEX(bb(w:), '"') - 2
                cc = bb(w:v)
                READ(cc,*) n
                nks = n
                EXIT
            ENDIF
        ENDDO
            
        !WRITE(0,*) 'Reading k-point Set Done'
          
        ! read nbnd
        n=0
        DO u = 1, 100
            READ(xml_pw_export,'(a)') bb
            if (INDEX(bb, 'nbnd="') /= 0) THEN
                w = INDEX(bb, 'nbnd="') + 6
                v = w + INDEX(bb(w:), '"') - 2
                cc = bb(w:v)
                READ(cc,*) n
                nbnd = n
                EXIT
            ENDIF
        ENDDO
    
        !WRITE(0,*) 'Reading Number of Band Done'
          
        ! read ngs
        n=0
        DO u = 1, 100
            READ(xml_pw_export,'(a)') bb
            IF (INDEX(bb, 'npw="') /= 0) THEN
                w = INDEX(bb, 'npw="') + 5
                v = w + INDEX(bb(w:), '"') - 2
                cc = bb(w:v)
                READ(cc,*) n
                ngs = n
                EXIT
            ENDIF
        ENDDO
    
        !WRITE(0,*) 'Reading Number of Main Grid Done'
          
        ! read alat
        r=0.0d0
        DO u = 1, 100
            READ(xml_pw_export,'(a)') bb
            IF (INDEX(bb, 'Data alat="') /= 0) THEN
                w = INDEX(bb, 'Data alat="') + 11
                v = w + INDEX(bb(w:), '"') - 2
                cc = bb(w:v)
                READ(cc,*) r
                alat = r
                EXIT
            ENDIF
        ENDDO
    
        !WRITE(0,*) 'Reading Alat Done'
          
        ! Read avec
        ALLOCATE(avec(3,3))
        DO u  = 1, 100
            READ(xml_pw_export,'(a)') bb
            IF (INDEX(bb, 'a1 xyz="') /= 0) THEN
                w = INDEX(bb, 'a1 xyz="') + 8
                v = w + INDEX(bb(w:), ' ') - 2
                cc = bb(w:v)
                READ(cc,*) r
                avec(1, 1) = r
          
                w = v + 2
                v = w + INDEX(bb(w:), ' ') - 2
                cc = bb(w:v)
                READ(cc,*) r
                avec(1, 2) = r
          
                w = v + 2
                v = w + INDEX(bb(w:), '"') - 2
                cc = bb(w:v)
                READ(cc,*) r
                avec(1, 3) = r
    
                !WRITE(0,*) 'Reading Primitive Lattice Vector a1 Done'
          
                EXIT
            ENDIF
        ENDDO
          
        DO u  = 1, 100
            READ(xml_pw_export,'(a)') bb
            IF (INDEX(bb, 'a2 xyz="') /= 0) THEN
                w = INDEX(bb, 'a2 xyz="') + 8
                v = w + INDEX(bb(w:), ' ') - 2
                cc = bb(w:v)
                READ(cc,*) r
                avec(2, 1) = r
          
                w = v + 2
                v = w + INDEX(bb(w:), ' ') - 2
                cc = bb(w:v)
                READ(cc,*) r
                avec(2, 2) = r
          
                w = v + 2
                v = w + INDEX(bb(w:), '"') - 2
                cc = bb(w:v)
                READ(cc,*) r
                avec(2, 3) = r
         
                !WRITE(0,*) 'Reading Primitive Lattice Vector a2 Done'
    
                EXIT
            ENDIF
        ENDDO
          
        DO u  = 1, 100
            READ(xml_pw_export,'(a)') bb
            IF (INDEX(bb, 'a3 xyz="') /= 0) THEN
                w = INDEX(bb, 'a3 xyz="') + 8
                v = w + INDEX(bb(w:), ' ') - 2
                cc = bb(w:v)
                READ(cc,*) r
                avec(3, 1) = r
          
                w = v + 2
                v = w + INDEX(bb(w:), ' ') - 2
                cc = bb(w:v)
                READ(cc,*) r
                avec(3, 2) = r
          
                w = v + 2
                v = w + INDEX(bb(w:), '"') - 2
                cc = bb(w:v)
                READ(cc,*) r
                avec(3, 3) = r
          
                !WRITE(0,*) 'Reading Primitive Lattice Vector a3 Done'
    
                EXIT
            ENDIF
        ENDDO
          
        ! Read bvec
        ALLOCATE(bvec(3,3))
        DO u  = 1, 100
            READ(xml_pw_export,'(a)') bb
            IF (INDEX(bb, 'b1 xyz="') /= 0) THEN
                w = INDEX(bb, 'b1 xyz="') + 8
                v = w + INDEX(bb(w:), ' ') - 2
                cc = bb(w:v)
                READ(cc,*) r
                bvec(1, 1) = r
          
                w = v + 2
                v = w + INDEX(bb(w:), ' ') - 2
                cc = bb(w:v)
                READ(cc,*) r
                bvec(1, 2) = r
          
                w = v + 2
                v = w + INDEX(bb(w:), '"') - 2
                cc = bb(w:v)
                READ(cc,*) r
                bvec(1, 3) = r
          
                !WRITE(0,*) 'Reading Reciprocal Lattice Vector b1 Done'
    
                EXIT
            ENDIF
        ENDDO
    
        DO u  = 1, 100
            READ(xml_pw_export,'(a)') bb
            IF (INDEX(bb, 'b2 xyz="') /= 0) THEN
                w = INDEX(bb, 'b2 xyz="') + 8
                v = w + INDEX(bb(w:), ' ') - 2
                cc = bb(w:v)
                READ(cc,*) r
                bvec(2, 1) = r
          
                w = v + 2
                v = w + INDEX(bb(w:), ' ') - 2
                cc = bb(w:v)
                READ(cc,*) r
                bvec(2, 2) = r
          
                w = v + 2
                v = w + INDEX(bb(w:), '"') - 2
                cc = bb(w:v)
                READ(cc,*) r
                bvec(2, 3) = r
          
                !WRITE(0,*) 'Reading Reciprocal Lattice Vector b2 Done'
    
                EXIT
            ENDIF
        ENDDO
    
        DO u  = 1, 100
            READ(xml_pw_export,'(a)') bb
            IF (INDEX(bb, 'b3 xyz="') /= 0) THEN
                w = INDEX(bb, 'b3 xyz="') + 8
                v = w + INDEX(bb(w:), ' ') - 2
                cc = bb(w:v)
                READ(cc,*) r
                bvec(3, 1) = r
                  
                w = v + 2
                v = w + INDEX(bb(w:), ' ') - 2
                cc = bb(w:v)
                READ(cc,*) r
                bvec(3, 2) = r
          
                w = v + 2
                v = w + INDEX(bb(w:), '"') - 2
                cc = bb(w:v)
                READ(cc,*) r
                bvec(3, 3) = r
    
                !WRITE(0,*) 'Reading Reciprocal Lattice Vector b3 Done'
          
                EXIT
            ENDIF
        ENDDO
          
        !READ weight
        ALLOCATE(wk(nks))
        sum_wk = 0.0d0
        DO
            READ(xml_pw_export,'(a)') bb
            IF (INDEX(bb, 'weights type="') /= 0) THEN
          
                DO i = 1, nks
                    READ(xml_pw_export, *) wk(i)
                    sum_wk = sum_wk + wk(i)
                ENDDO
          
                EXIT
            ENDIF
        ENDDO
        !WRITE(0,*) 'Reading Weight Done'
          
        !Read k
        ALLOCATE(k(nks, 3))
        DO
            READ(xml_pw_export,'(a)') bb
            IF (INDEX(bb, '<k type="') /= 0) THEN
          
                DO i = 1, nks
                    READ(xml_pw_export, *) k(i, 1), k(i, 2), k(i, 3)
                ENDDO
          
                EXIT
            ENDIF
        ENDDO
    
        WRITE(0,*) 'Reading k-point Done'
          
        k(:,:) = k(:,:)*(2.0d0*pi/alat)
          
        !READ G
        ALLOCATE(g2(ngs, 3))
        ALLOCATE(gi(ngs, 3))
        DO
            READ(xml_pw_export, '(a)') bb
            IF (INDEX(bb, '<g type="')) THEN
                DO i = 1, ngs
                    READ(xml_pw_export, *) gi(i, 1), gi(i, 2), gi(i, 3)
                    !write(0,*) g1, g2, g3
                    g2(i, 1) = dble(gi(i, 1))*bvec(1,1) + dble(gi(i, 2))*bvec(2,1) + dble(gi(i, 3))*bvec(3,1)
                    g2(i, 2) = dble(gi(i, 1))*bvec(1,2) + dble(gi(i, 2))*bvec(2,2) + dble(gi(i, 3))*bvec(3,2)
                    g2(i, 3) = dble(gi(i, 1))*bvec(1,3) + dble(gi(i, 2))*bvec(2,3) + dble(gi(i, 3))*bvec(3,3)
                    !WRITE(1025,*) gi(i, 1), gi(i, 2), gi(i, 3)
                ENDDO
                EXIT
            ENDIF
        ENDDO
    
        WRITE(0,*) 'Reading G Done'
          
        !Read G for each K-point
        allocate(gik(nks,ngs), n_gik(nks))
        gik(:,:) = 0
        DO
            READ(xml_pw_export, '(a)') bb
            IF (INDEX(bb, '<Wfc_grids')) THEN
                i = 1
                DO
                    READ(xml_pw_export, '(a)') bb
                    IF (INDEX(bb, '<Kpoint.')) THEN
                        READ(xml_pw_export, '(a)') bb
                        w = INDEX(bb, 'size="') + 6
                        v = w + INDEX(bb(w:), '"') - 2
                        cc = bb(w:v)
                        READ(cc,*) n_gik(i)
          
                        DO u = 1, n_gik(i)
                            READ(xml_pw_export, *) gik(i,u)
                        ENDDO
          
                        i = i + 1
                    ENDIF

                    IF (i > nks) THEN
                        EXIT
                    ENDIF
                ENDDO
                EXIT
            ENDIF
        ENDDO
    
        !WRITE(0,*) 'Reading G for k-points Done'
          
        ALLOCATE(g(1:maxval(gik),3))
        g(:,:) = g2(1:maxval(n_gik),:)
        deallocate(g2)
        ngs = maxval(gik)
          
        !Read eigv
        allocate(eigv(nks, nbnd))
        DO
            READ(xml_pw_export, '(a)') bb
            IF (INDEX(bb, '<Eigenvalues ')) THEN
                DO u = 1, nks
                    READ(xml_pw_export, '(a)') bb
          
                    DO v = 1, nbnd
                        READ(xml_pw_export, *) eigv(u, v)                   
                    ENDDO
          
                    READ(xml_pw_export, '(a)') bb
                ENDDO
                EXIT
            ENDIF
        ENDDO
    
        WRITE(0,*) 'Reading Eigenvalues Done'
          
        !WRITE(0,*) 'Starting Allocate Eigenvector'
        !eigv(:,:) = eigv(:,:)*ry2ev                                    ! Rydberg --> eV
        allocate(coef_wfc(nks, nbnd, ngs, 1:nspin))
        !WRITE(0,*) 'Allocate Eigenvector Done'
    
        !WRITE(0,*) 'Starting Initial Zero Matrix'
        coef_wfc(:,:,:,:) = 0.0d0
        !WRITE(0,*) 'Initial Zero Matrix Done'
          
        !Read c
        !WRITE(0,*) 'Starting Eigenvector'
        DO
            READ(xml_pw_export, '(a)') bb
            IF (INDEX(bb, '<Eigenvectors')) THEN
                DO i = 1, nks
                    DO 
                    READ(xml_pw_export, '(a)') bb
                        IF (INDEX(bb, '<Wfc.')) THEN
                            DO v = 1, nbnd
                            
                                DO u = 1, n_gik(i)
                                    READ(xml_pw_export,*) r, s
                                    coef_wfc(i, v, gik(i,u),ispin) = cmplx(r,s)
                                ENDDO
          
                            READ(xml_pw_export, '(a)') bb
                            READ(xml_pw_export, '(a)') bb
    
                            ENDDO
                            EXIT
                        ENDIF
                    ENDDO
                ENDDO
                EXIT
            ENDIF
        ENDDO
    
        WRITE(0,*) 'Reading Eigenvector Done'
    
        close(xml_pw_export)
    
        WRITE(0,*) 'Reading XML File Program Done'
            
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE reading_xmlfile
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    
!----------------------------------------------------------------------------------------------------------------------------------------------------------
END MODULE readingfile_mod
!----------------------------------------------------------------------------------------------------------------------------------------------------------
    
!----------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE dvec_mod
!----------------------------------------------------------------------------------------------------------------------------------------------------------
    
    USE inputoutput_mod
    USE dipolevector_mod , ONLY : prefix, outdir, input_pw_export, &
                                sorb, valence_band, conduction_band, norm, &
                                nks, nbnd, ngs, nspin, npol, &
                                k, g, wk, eigv, &
                                dvec, polvec, coef_wfc, &
                                matele_opt, &
                                pi, im
    
    IMPLICIT NONE
    
    CONTAINS
       
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE electronphoton_interaction
    !------------------------------------------------------------------------------------------------------------------------------------------------------
       
    ! This subroutine calculate the electron photon matrix element
    INTEGER :: ik, ibi, ibf, ig, ip, ispin, ixyz
    INTEGER :: nev, nec
    ! Band index for dipole vector without considering spin, or with considering spin by taking the bottom one 
    INTEGER :: nev2, nec2
    ! Band index for dipole vector with considering spin by taking the upper one

    REAL(KIND=KIND(1.0d0)) :: eps=1.0d-6
    
    COMPLEX (KIND=KIND(1.0d0)) :: coeff, phase                          ! INTEGER :: i_max_g
    COMPLEX (KIND=KIND(1.0d0)) :: max_g                                 ! REAL(KIND=KIND(1.0d0)), ALLOCATABLE :: matele_opt2(:,:,:,:) 
    
    nev = valence_band                                                  ! Valence Band
    nec = conduction_band                                               ! Conduction Band

    !WRITE(0,*) 'Valence Band Index :', nev
    !WRITE(0,*) 'Conduction Band Index :', nec

    if (sorb .eq. .true.) then
        nev2 = nev + 1
        nec2 = nec + 1

        !WRITE(0,*) 'Considering Spin Orbit Coupling'
        !WRITE(0,*) 'Valence Band Index 2 :', nev2
        !WRITE(0,*) 'Conduction Band Index 2 :', nec2
    end if

    ALLOCATE(dvec(nks,nbnd,nbnd,3))                                     
    dvec(:,:,:,:) = (0.0d0,0.0d0)
    
    WRITE(0,*) 'Calling Dipole Vector Formula for Calculation'
    
    DO ik = 1, nks
        DO ibi = 1, nbnd                                                ! Index of energy bands of initial state
            DO ibf = 1, nbnd                                            ! Index of energy bands of final state

                ! Calculation of dipole vector
                DO ispin = 1, nspin
                    DO ig = 1, ngs
                    coeff = im * coef_wfc(ik,ibi,ig,ispin) * CONJG(coef_wfc(ik,ibf,ig,ispin))
                        DO ixyz = 1, 3
                            dvec(ik,ibi,ibf,ixyz) = coeff * (k(ik,ixyz) + g(ig,ixyz)) + dvec(ik,ibi,ibf,ixyz)                     
                        END DO
                    END DO
                END DO

            END DO
        END DO
    END DO

    !WRITE(0,*) 'Dipole Vector Formula Done'
    WRITE(0,*) 'Starting Dipole Vector Calculation'

    ! xml_read
    do ik = 1, nks
        do ibi = 1, nbnd
            do ibf = 1, nbnd
                dvec(ik,ibi,ibf,:) = dvec(ik,ibi,ibf,:)*(cos(-atan2(aimag(dvec(ik,ibi,ibf,2)),real(dvec(ik,ibi,ibf,2)))) + im*sin(-atan2(aimag(dvec(ik,ibi,ibf,2)),real(dvec(ik,ibi,ibf,2)))))
            end do
        end do
    end do
    
    !WRITE(0,*) 'Dipole Vector Calculation Done'
    WRITE(0,*) 'Writting Kpoint Data'
    WRITE(0,*) 'Writting Dipole Vector and Oscillator Strength Data'

    ! For xml_reading dipole vector
    open(86, file = "kpoint.dat")
    open(87, file = "dipole_vector.dat")                                ! Dipole vector without considering spin, or with considering spin by taking the bottom one
    open(88, file = "oscillator_strength.dat")                          ! Oscillator strength without considering spin, or with considering spin by taking the bottom one

    write(86,*) '# k(x) k(y) k(z)'
    write(86,*) '# If you need this file for another program input, delete all of this messages'

    write(87,*) '# valence =', nev, ', conduction =', nec
    write(87,*) '# k(x) k(y) k(z) Real(x) Real(y) Real(z) Imaginary(x) Imaginary(y) Imaginary(z)'
    write(87,*) '# If you need this file for another program input, delete all of this messages'

    write(88,*) '# valence =', nev, ', conduction =', nec
    write(88,*) '# k(x) k(y) k(z) Oscillator Strength'
    write(88,*) '# If you need this file for another program input, delete all of this messages'

    do ik = 1, nks
        write(86,'(3(1x,f10.4))') k(ik,1), k(ik,2), k(ik,3)
        write(87,'(3(1x,f10.4),6(1x,e16.8e3))') k(ik,1), k(ik,2), k(ik,3), real(dvec(ik,nev,nec,1))*norm, real(dvec(ik,nev,nec,2))*norm, real(dvec(ik,nev,nec,3))*norm, aimag(dvec(ik,nev,nec,1))*norm, aimag(dvec(ik,nev,nec,2))*norm, aimag(dvec(ik,nev,nec,3))*norm
        write(88,'(3(1x,f10.4),1(1x,e16.8e3))') k(ik,1), k(ik,2), k(ik,3), sqrt(abs(real(dvec(ik,nev,nec,1))*aimag(dvec(ik,nev,nec,1)) + real(dvec(ik,nev,nec,2))*aimag(dvec(ik,nev,nec,2))))
    end do

    write(86,*) ' '
    write(87,*) ' '
    write(88,*) ' '

    close(86)
    close(87)
    close(88)

    if (sorb .eq. .true.) then
        open(95, file = "dipole_vector2.dat")                           ! Dipole vector with considering spin by taking the upper one
        open(96, file = "oscillator_strength2.dat")                     ! Oscillator strength with considering spin by taking the upper one

        write(95,*) '# valence =', nev2, ', conduction =', nec2
        write(95,*) '# k(x) k(y) k(z) Real(x) Real(y) Real(z) Imaginary(x) Imaginary(y) Imaginary(z)'
        write(95,*) '# If you need this file for another program input, delete all of this messages'

        write(96,*) '# valence =', nev2, ', conduction =', nec2
        write(96,*) '# k(x) k(y) k(z) Oscillator Strength'
        write(96,*) '# If you need this file for another program input, delete all of this messages'

        do ik = 1, nks
            write(95,'(3(1x,f10.4),6(1x,e16.8e3))') k(ik,1), k(ik,2), k(ik,3), real(dvec(ik,nev2,nec2,1))*norm, real(dvec(ik,nev2,nec2,2))*norm, real(dvec(ik,nev2,nec2,3))*norm, aimag(dvec(ik,nev2,nec2,1))*norm, aimag(dvec(ik,nev2,nec2,2))*norm, aimag(dvec(ik,nev2,nec2,3))*norm
            write(96,'(3(1x,f10.4),1(1x,e16.8e3))') k(ik,1), k(ik,2), k(ik,3), sqrt(abs(real(dvec(ik,nev2,nec2,1))*aimag(dvec(ik,nev2,nec2,1)) + real(dvec(ik,nev2,nec2,2))*aimag(dvec(ik,nev2,nec2,2))))
        end do

        write(95,*) ' '
        write(96,*) ' '

        close(95)
        close(96)
    end if

    !stop('finish')

    !WRITE(0,*) 'Writting Kpoint Data Done'
    !WRITE(0,*) 'Writting Dipole Vector and Oscillator Strength Data Done'
    
    DEALLOCATE(coef_wfc, g)

    WRITE(0,*) 'Dipole Vector Program Done'
    
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE electronphoton_interaction
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    
!----------------------------------------------------------------------------------------------------------------------------------------------------------
END MODULE dvec_mod
!----------------------------------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE matele_mod
!----------------------------------------------------------------------------------------------------------------------------------------------------------
    
    USE inputoutput_mod
    USE dipolevector_mod
    USE readingfile_mod

    IMPLICIT NONE

    CONTAINS

    !------------------------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE out_file
    !------------------------------------------------------------------------------------------------------------------------------------------------------
            
    INTEGER :: ik
    
    WRITE(stdout, '(1x)')
    
    WRITE(stdout, '(5x,"Number of k points",I6)') nks
    WRITE(stdout, '(5x,"Number of bands",I6)') nbnd
    WRITE(stdout, '(5x,"Number of G vectors",I6)') ngs
        
    WRITE(stdout, '(1x)')
    WRITE(stdout, '(5x,"Lattice constant alat = ",f8.4,3x,"(au)")') alat
    WRITE(stdout, '(1x)')
    
    WRITE(stdout,'(5x,"Primitive lattice vectors (au)")')
    WRITE(stdout,'(6x,A6,3(1x,f8.4),A1)') 'a1 = (', avec(1,1), avec(1,2), avec(1,3), ')'
    WRITE(stdout,'(6x,A6,3(1x,f8.4),A1)') 'a2 = (', avec(2,1), avec(2,2), avec(2,3), ')'
    WRITE(stdout,'(6x,A6,3(1x,f8.4),A1)') 'a3 = (', avec(3,1), avec(3,2), avec(3,3), ')'
    WRITE(stdout,'(5x,"Reciprocal lattice vectors (au)")')
    WRITE(stdout,'(6x,A6,3(1x,f8.4),A1)') 'b1 = (', bvec(1,1), bvec(1,2), bvec(1,3), ')'
    WRITE(stdout,'(6x,A6,3(1x,f8.4),A1)') 'b2 = (', bvec(2,1), bvec(2,2), bvec(2,3), ')'
    WRITE(stdout,'(6x,A6,3(1x,f8.4),A1)') 'b3 = (', bvec(3,1), bvec(3,2), bvec(3,3), ')'
    
    WRITE(stdout, '(1x)')
    
    WRITE(stdout,'(5x,"k points (au)")')
    DO ik = 1, nks
        WRITE(stdout,'(6x,I5,3x,A1,3(1x,f8.4),A1,1x,f8.4)') ik,'(', k(ik,1), k(ik,2), k(ik,3), ')', wk(ik)
    ENDDO
    
    WRITE(stdout, '(1x)')
    
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE out_file
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE matrix_element_transition
    !------------------------------------------------------------------------------------------------------------------------------------------------------

    ! Output for data files
    INTEGER :: ik, ibi, ibf, ip, ibnd, inb

    REAL :: nb
    REAL(KIND=KIND(1.0d0)) :: eps=1.0d-6
    
    CHARACTER(len=256) :: make_outdir, make_outdir_matele_opt, make_outdir_eigenvalues
    CHARACTER(len=256) :: outdir_matele_opt, outdir_eigenvalues
    CHARACTER(len=256) :: filename_matele_opt, readme_mopt
    CHARACTER(len=256) :: filename_eigenvalues, readme_eigv
    CHARACTER(len=256) :: ecv, cmopt1, cmopt2                            ! Band index transition

    nb = nbnd*nbnd

    WRITE(make_outdir,'("mkdir -p ",A228)') outdir
    CALL SYSTEM(make_outdir)

    ! Data output for eigen energy list
    outdir_eigenvalues = TRIM(ADJUSTL(outdir))//'/eigenvalues/'
    WRITE(make_outdir_eigenvalues,'("mkdir -p ",A228)') outdir_eigenvalues
    CALL SYSTEM(make_outdir_eigenvalues)

    readme_eigv = TRIM(ADJUSTL(outdir_eigenvalues))//'00-readme.txt'
    OPEN(32, file=readme_eigv)

    WRITE(32,*) 'Eigenvalue data have been written'
    WRITE(32,*) ' '
    WRITE(32,*) 'Filename are eigenvalues_(band).dat'
    WRITE(32,*) 'There are ', nbnd, ' data files due to eigenvalues correspond to each band'
    WRITE(32,*) 'Eigenvalue data structure : '
    write(32,*) 'k(x) k(y) k(z) Eigenvalues'
    WRITE(32,*) ' '
    WRITE(32,*) 'End of structure'

    WRITE(0,*) 'Make Eigenvalues Directory Done'
    WRITE(0,*) 'Writting Eigenvalues Data'

    DO inb = 1, nbnd
        WRITE(ecv,'(I3)') inb
        filename_eigenvalues = TRIM(ADJUSTL(outdir_eigenvalues))//'eigenvalues_'//TRIM(ADJUSTL(ecv))//'.dat'
        OPEN(fn_eigv, file=filename_eigenvalues)

        DO ik = 1, nks
            WRITE(fn_eigv,'(3(1x,f10.4),1(1x,e16.8e3))') k(ik,1), k(ik,2), k(ik,3), eigv(ik,inb)*ry2ev
        ENDDO
    ENDDO

    write(fn_eigv,*) ' '

    !WRITE(0,*) 'Writting Eigenvalues Data Done'

    ! Data output for electron-photon matrix element
    IF (plot_matele_opt .EQ. .TRUE.) THEN
        WRITE(0,*) 'Calling Matrix Element Transition for Calculation'

        ALLOCATE(matele_opt(nks,nbnd,nbnd,npol))
        matele_opt(:,:,:,:) = 0.0d0
        outdir_matele_opt = TRIM(ADJUSTL(outdir))//'/matele_opt/'
        WRITE(make_outdir_matele_opt,'("mkdir -p ",A228)') outdir_matele_opt
        CALL SYSTEM(make_outdir_matele_opt)

        readme_mopt = TRIM(ADJUSTL(outdir_matele_opt))//'00-readme.txt'
        OPEN(35, file=readme_mopt)

        WRITE(35,*) 'Optical matrix element data have been written'
        WRITE(35,*) ' '
        WRITE(35,*) 'Filename are matele_opt_(ibi)_(ibf).dat'
        WRITE(35,*) 'With ibi is the index of initial band, ibf is index of final band'
        WRITE(35,*) 'There are', nb, 'data files due to excitation from ibi to ibf band'
        WRITE(35,*) 'Optical matrix element data structure : '

        IF (circular_pol .EQ. .FALSE. .AND.  nonpol .EQ. .FALSE.) THEN
            write(35,*) 'k(x) k(y) k(z) Absolute(x) Absolute(y) Absolute(z)'
        ELSE IF (circular_pol .EQ. .TRUE. .AND.  nonpol .EQ. .FALSE.) THEN
            write(35,*) 'k(x) k(y) k(z) Absolute(Left)^2 Absolute(Right)^2'
        ELSE                           
            write(35,*) 'k(x) k(y) k(z) Absolute(non)'
        END IF

        WRITE(35,*) ' '
        WRITE(35,*) 'End of structure'

        WRITE(0,*) 'Make Electron Photon Matrix Element Transition Directory Done'
        WRITE(0,*) 'Starting Electron Photon Matrix Element Transition Calculation'
        WRITE(0,*) 'Writting Electron Photon Matrix Element Transition Data'

        DO ibi = 1, nbnd
            WRITE(cmopt1,'(I3)') ibi

            DO ibf = 1, nbnd
                WRITE(cmopt2,'(I3)') ibf
                filename_matele_opt = TRIM(ADJUSTL(outdir_matele_opt))//'matele_opt_'//TRIM(ADJUSTL(cmopt1))//'_'//TRIM(ADJUSTL(cmopt2))//'.dat'
                OPEN(fn_matele_opt, file=filename_matele_opt)

                DO ik = 1, nks

                    DO ip = 1, npol
                        matele_opt(ik,ibi,ibf,ip) = polvec(ip,1)*CONJG(dvec(ik,ibi,ibf,1)) + polvec(ip,2)*CONJG(dvec(ik,ibi,ibf,2)) + polvec(ip,3)*CONJG(dvec(ik,ibi,ibf,3))
                    END DO
                    
                    IF (circular_pol .EQ. .FALSE. .AND.  nonpol .EQ. .FALSE.) THEN
                        WRITE(fn_matele_opt,'(3(1x,f10.4),9(1x,e16.8e3))') k(ik,1), k(ik,2), k(ik,3), (ABS(matele_opt(ik,ibi,ibf,1)))**2, ABS(matele_opt(ik,ibi,ibf,2))**2, ABS(matele_opt(ik,ibi,ibf,3))**2
                    ELSE IF (circular_pol .EQ. .TRUE. .AND.  nonpol .EQ. .FALSE.) THEN
                        WRITE(fn_matele_opt,'(3(1x,f10.4),6(1x,e16.8e3))') k(ik,1), k(ik,2), k(ik,3), ABS(matele_opt(ik,ibi,ibf,4))**2, ABS(matele_opt(ik,ibi,ibf,5))**2
                    ELSE                           
                        WRITE(fn_matele_opt,'(3(1x,f10.4),15(1x,e16.8e3))') k(ik,1), k(ik,2), k(ik,3), ABS(matele_opt(ik,ibi,ibf,6))**2
                    END IF

                END DO
            END DO
        END DO

        write(fn_matele_opt,*) ' '

        !WRITE(0,*) 'Writting Electron Photon Matrix Element Transition Done'
        !WRITE(0,*) 'Electron Photon Matrix Element Transition Done'
    END IF

    WRITE(0,*) 'Matrix Element Transition Program Done'

    !------------------------------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE matrix_element_transition
    !------------------------------------------------------------------------------------------------------------------------------------------------------
    
!----------------------------------------------------------------------------------------------------------------------------------------------------------
END MODULE matele_mod
!----------------------------------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------------------------------
PROGRAM mainprogram
!----------------------------------------------------------------------------------------------------------------------------------------------------------

    USE inputoutput_mod, ONLY : stdout
    USE dipolevector_mod
    USE readingfile_mod
    USE dvec_mod
    USE matele_mod

    IMPLICIT NONE
    
    ! For date and time
    integer :: ti(1:8), i
    character (len=10) :: t(1:3)
    
    ! Initial message
    WRITE(stdout,'(/5x,a)') repeat('=',75)
    WRITE(stdout,'(a)')"     "
    WRITE(stdout,'(a)')"     ----------------------------------------------------------------------"
    WRITE(stdout,'(a)')"                Main Program Optical Matrix Element Transition"
    WRITE(stdout,'(a)')"     ----------------------------------------------------------------------"
    WRITE(stdout,'(a)')"     "
    WRITE(stdout,'(a)')"     Main Program  for Optical Matrix Element Transition"
    WRITE(stdout,'(a)')"     Created by Yuki Tatsumi, started from 2016/11/15"
    WRITE(stdout,'(a)')"     Update by R. Layton & MJ. Prakoso, started from 2020/09/04"
    WRITE(stdout,'(a)')"     "
    WRITE(stdout,'(/5x,a)') repeat('=',75)
    
    CALL date_and_time(t(1),t(2),t(3),ti)
    WRITE(stdout,'(A,i0,A,i0,A,i0,A,i0,A,i0,A,i0)') "     Program started at ", &
    ti(1),"/",ti(2),"/",ti(3)," , ",ti(5),":",ti(6),":",ti(7)
    !WRITE(0,*) 'Using Date and Time for Initial State Done'

    WRITE(0,'(A,i0,A,i0,A,i0,A,i0,A,i0,A,i0)') " Program started at ", &
    ti(1),"/",ti(2),"/",ti(3)," , ",ti(5),":",ti(6),":",ti(7)

    ! Initial setting for parameters
    CALL unit_parameter
    !WRITE(0,*) 'Using Parameter Done'
    
    ! Reading the input file
    CALL reading_QEinputfile
    !WRITE(0,*) 'Using Reading Input Program Done'
    
    ! Reading the data file of wave function etc. calculated by Quantum Espresso
    CALL reading_xmlfile
    !WRITE(0,*) 'Using Reading XML File Program Done'
    
    ! Initial setting of polarization vectors
    CALL polarization_vector
    !WRITE(0,*) 'Using Polarized Vector Program Done'

    ! Calculation and obtaining the dipole vectors
    CALL electronphoton_interaction
    !WRITE(0,*) 'Using Dipole Vector Program Done'

    ! Standard output of reading data
    CALL out_file
    !WRITE(0,*) 'Using Standart Output File Program Done'

    ! Calculation and obtaining the electron photon matrix element transition
    CALL matrix_element_transition
    !WRITE(0,*) 'Using Matrix Element Transition Program Done'
    
    CALL date_and_time(t(1),t(2),t(3),ti)
    WRITE(stdout, '(1x)')
    WRITE(stdout,'(A,i0,A,i0,A,i0,A,i0,A,i0,A,i0)') "     Program finished at ", &
    ti(1),"/",ti(2),"/",ti(3)," , ",ti(5),":",ti(6),":",ti(7) 

    !WRITE(0,*) 'Using Date and Time for Final State Done'
    
    WRITE(stdout,'(/5x,a)') repeat('=',75)
    WRITE(stdout,'(5x,"JOB DONE")')
    WRITE(stdout,'(/5x,a)') repeat('=',75)

    WRITE(0,*) 'Job Done'
    WRITE(0,'(A,i0,A,i0,A,i0,A,i0,A,i0,A,i0)') " Program finished at ", &
    ti(1),"/",ti(2),"/",ti(3)," , ",ti(5),":",ti(6),":",ti(7)
    WRITE(0,*) ''
    
!----------------------------------------------------------------------------------------------------------------------------------------------------------
END PROGRAM mainprogram
!----------------------------------------------------------------------------------------------------------------------------------------------------------