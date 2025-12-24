!==============================================================================
! Emulator Input/Output Module
! Handles parsing of EMULATOR_TRAIN and EMULATOR_PREDICT namelists
!==============================================================================
module emulator_io

use precision
implicit none

integer, parameter :: MAX_PARAMS = 20
integer, parameter :: MAX_PARAM_SETS = 1000

!------------------------------------------------------------------------------
! Parameter range type (for training)
!------------------------------------------------------------------------------
type :: param_range_t
    character(len=32) :: name = ''      ! e.g., 't_1_uv'
    real(dpreal) :: min_val = 0.0_dpreal
    real(dpreal) :: max_val = 0.0_dpreal
    ! Parsed components
    character(len=4) :: kp1 = ''        ! Potential type: 't', 'x', 'b', 'p', 'a'
    integer :: kp2 = 0                   ! Potential index
    character(len=8) :: param = ''       ! Parameter: 'uv', 'wd', etc.
end type

!------------------------------------------------------------------------------
! Parameter set type (for prediction)
!------------------------------------------------------------------------------
type :: param_set_t
    integer :: n_params = 0
    character(len=32) :: names(MAX_PARAMS)
    real(dpreal) :: values(MAX_PARAMS)
end type

!------------------------------------------------------------------------------
! Training configuration
!------------------------------------------------------------------------------
type :: train_config_t
    character(len=256) :: base_input = ''
    integer :: n_samples = 50
    integer :: n_basis = 0              ! 0 = auto (GS method)
    character(len=8) :: method = 'gs'   ! 'svd' or 'gs'
    real(dpreal) :: tol = 1.0d-8
    character(len=256) :: output_file = 'emulator.dat'
    logical :: multi_J = .true.         ! Train over all J values (J as parameter)
    logical :: random_J = .false.       ! Randomly sample J (vs all J per sample)
    integer :: n_J_samples = 0          ! Number of J samples per parameter set (0=all)
    ! Parameter ranges
    integer :: n_params = 0
    type(param_range_t) :: params(MAX_PARAMS)
end type

!------------------------------------------------------------------------------
! Prediction configuration
!------------------------------------------------------------------------------
type :: predict_config_t
    character(len=256) :: emulator_file = ''
    character(len=256) :: base_input = ''
    character(len=16) :: predict_method = 'matrix'  ! 'rbf' or 'matrix'
    character(len=32) :: output_format = 'cross_section'
    ! Parameter sets
    integer :: n_sets = 0
    type(param_set_t) :: sets(MAX_PARAM_SETS)
end type

contains

!==============================================================================
! Parse parameter name (e.g., 't_1_uv' -> kp1='t', kp2=1, param='uv')
!==============================================================================
subroutine parse_param_name(name, kp1, kp2, param, ierr)
    character(len=*), intent(in) :: name
    character(len=*), intent(out) :: kp1, param
    integer, intent(out) :: kp2, ierr

    integer :: i1, i2
    character(len=32) :: temp

    ierr = 0
    kp1 = ''
    kp2 = 0
    param = ''

    ! Find first underscore
    i1 = index(name, '_')
    if (i1 == 0) then
        ierr = 1
        return
    end if

    ! Extract kp1
    kp1 = name(1:i1-1)

    ! Find second underscore
    temp = name(i1+1:)
    i2 = index(temp, '_')
    if (i2 == 0) then
        ierr = 2
        return
    end if

    ! Extract kp2
    read(temp(1:i2-1), *, iostat=ierr) kp2
    if (ierr /= 0) return

    ! Extract param name
    param = trim(temp(i2+1:))

end subroutine

!==============================================================================
! Read training configuration
!==============================================================================
subroutine read_train_config(filename, config, ierr)
    character(len=*), intent(in) :: filename
    type(train_config_t), intent(out) :: config
    integer, intent(out) :: ierr

    ! Namelist variables
    character(len=256) :: base_input, output_file
    integer :: n_samples, n_basis, n_J_samples
    character(len=8) :: method
    real(dpreal) :: tol
    logical :: multi_J, random_J

    character(len=32) :: name
    real(dpreal) :: min_val, max_val

    integer :: iunit, ios, ip

    namelist /EMULATOR_TRAIN/ base_input, n_samples, n_basis, method, tol, output_file, &
                              multi_J, random_J, n_J_samples
    namelist /PARAM_RANGE/ name, min_val, max_val

    ierr = 0
    iunit = 99

    ! Initialize defaults
    base_input = ''
    n_samples = 50
    n_basis = 0
    method = 'gs'
    tol = 1.0d-8
    output_file = 'emulator.dat'
    multi_J = .true.
    random_J = .false.
    n_J_samples = 0

    ! Open file
    open(unit=iunit, file=trim(filename), status='old', iostat=ios)
    if (ios /= 0) then
        ierr = -1
        write(*,'(A)') 'ERROR: Cannot open file: '//trim(filename)
        return
    end if

    ! Read EMULATOR_TRAIN namelist
    read(iunit, nml=EMULATOR_TRAIN, iostat=ios)
    if (ios /= 0) then
        ierr = -2
        write(*,'(A)') 'ERROR: Failed to read EMULATOR_TRAIN namelist'
        close(iunit)
        return
    end if

    config%base_input = base_input
    config%n_samples = n_samples
    config%n_basis = n_basis
    config%method = method
    config%tol = tol
    config%output_file = output_file
    config%multi_J = multi_J
    config%random_J = random_J
    config%n_J_samples = n_J_samples

    ! Read PARAM_RANGE namelists
    config%n_params = 0
    do
        name = ''
        min_val = 0.0_dpreal
        max_val = 0.0_dpreal

        read(iunit, nml=PARAM_RANGE, iostat=ios)
        if (ios /= 0) exit
        if (len_trim(name) == 0) exit

        config%n_params = config%n_params + 1
        ip = config%n_params

        if (ip > MAX_PARAMS) then
            ierr = -3
            write(*,'(A,I3)') 'ERROR: Too many parameters, max=', MAX_PARAMS
            exit
        end if

        config%params(ip)%name = name
        config%params(ip)%min_val = min_val
        config%params(ip)%max_val = max_val

        ! Parse parameter name
        call parse_param_name(name, config%params(ip)%kp1, &
                             config%params(ip)%kp2, &
                             config%params(ip)%param, ios)
        if (ios /= 0) then
            write(*,'(A)') 'WARNING: Failed to parse parameter name: '//trim(name)
        end if
    end do

    close(iunit)

    ! Print configuration
    call print_train_config(config)

end subroutine

!==============================================================================
! Read training configuration from stdin
!==============================================================================
subroutine read_train_config_stdin(config, ierr)
    type(train_config_t), intent(out) :: config
    integer, intent(out) :: ierr

    ! Namelist variables
    character(len=256) :: base_input, output_file
    integer :: n_samples, n_basis, n_J_samples
    character(len=8) :: method
    real(dpreal) :: tol
    logical :: multi_J, random_J

    character(len=32) :: name
    real(dpreal) :: min_val, max_val

    integer :: ios, ip

    namelist /EMULATOR_TRAIN/ base_input, n_samples, n_basis, method, tol, output_file, &
                              multi_J, random_J, n_J_samples
    namelist /PARAM_RANGE/ name, min_val, max_val

    ierr = 0

    ! Initialize defaults
    base_input = ''
    n_samples = 50
    n_basis = 0
    method = 'gs'
    tol = 1.0d-8
    output_file = 'emulator.dat'
    multi_J = .true.
    random_J = .false.
    n_J_samples = 0

    ! Read EMULATOR_TRAIN namelist from stdin
    read(*, nml=EMULATOR_TRAIN, iostat=ios)
    if (ios /= 0) then
        ierr = -2
        write(*,'(A)') 'ERROR: Failed to read EMULATOR_TRAIN namelist from stdin'
        return
    end if

    config%base_input = base_input
    config%n_samples = n_samples
    config%n_basis = n_basis
    config%method = method
    config%tol = tol
    config%output_file = output_file
    config%multi_J = multi_J
    config%random_J = random_J
    config%n_J_samples = n_J_samples

    ! Read PARAM_RANGE namelists
    config%n_params = 0
    do
        name = ''
        min_val = 0.0_dpreal
        max_val = 0.0_dpreal

        read(*, nml=PARAM_RANGE, iostat=ios)
        if (ios /= 0) exit
        if (len_trim(name) == 0) exit

        config%n_params = config%n_params + 1
        ip = config%n_params

        if (ip > MAX_PARAMS) then
            ierr = -3
            write(*,'(A,I3)') 'ERROR: Too many parameters, max=', MAX_PARAMS
            exit
        end if

        config%params(ip)%name = name
        config%params(ip)%min_val = min_val
        config%params(ip)%max_val = max_val

        ! Parse parameter name
        call parse_param_name(name, config%params(ip)%kp1, &
                             config%params(ip)%kp2, &
                             config%params(ip)%param, ios)
        if (ios /= 0) then
            write(*,'(A)') 'WARNING: Failed to parse parameter name: '//trim(name)
        end if
    end do

    ! Print configuration
    call print_train_config(config)

end subroutine

!==============================================================================
! Read prediction configuration
!==============================================================================
subroutine read_predict_config(filename, config, ierr)
    character(len=*), intent(in) :: filename
    type(predict_config_t), intent(out) :: config
    integer, intent(out) :: ierr

    ! Namelist variables
    character(len=256) :: emulator_file, base_input
    character(len=16) :: predict_method
    character(len=32) :: output_format

    ! PARAM_SET variables - all 18 optical potential parameters
    real(dpreal) :: t_1_uv, t_1_rv, t_1_av  ! proton volume real
    real(dpreal) :: t_1_uw, t_1_rw, t_1_aw  ! proton volume imag
    real(dpreal) :: t_1_wd, t_1_rwd, t_1_awd ! proton surface imag
    real(dpreal) :: x_1_uv, x_1_rv, x_1_av  ! neutron volume real
    real(dpreal) :: x_1_uw, x_1_rw, x_1_aw  ! neutron volume imag
    real(dpreal) :: x_1_wd, x_1_rwd, x_1_awd ! neutron surface imag
    real(dpreal) :: b_1_uv, b_1_wd, p_1_uv  ! other potentials

    integer :: iunit, ios, is, ip
    logical :: has_value

    namelist /EMULATOR_PREDICT/ emulator_file, base_input, predict_method, output_format
    namelist /PARAM_SET/ t_1_uv, t_1_rv, t_1_av, t_1_uw, t_1_rw, t_1_aw, t_1_wd, t_1_rwd, t_1_awd, &
                         x_1_uv, x_1_rv, x_1_av, x_1_uw, x_1_rw, x_1_aw, x_1_wd, x_1_rwd, x_1_awd, &
                         b_1_uv, b_1_wd, p_1_uv

    ierr = 0
    iunit = 99

    ! Initialize defaults
    emulator_file = ''
    base_input = ''
    predict_method = 'matrix'
    output_format = 'cross_section'

    ! Open file
    open(unit=iunit, file=trim(filename), status='old', iostat=ios)
    if (ios /= 0) then
        ierr = -1
        write(*,'(A)') 'ERROR: Cannot open file: '//trim(filename)
        return
    end if

    ! Read EMULATOR_PREDICT namelist
    read(iunit, nml=EMULATOR_PREDICT, iostat=ios)
    if (ios /= 0) then
        ierr = -2
        write(*,'(A)') 'ERROR: Failed to read EMULATOR_PREDICT namelist'
        close(iunit)
        return
    end if

    config%emulator_file = emulator_file
    config%base_input = base_input
    config%predict_method = predict_method
    config%output_format = output_format

    ! Read PARAM_SET namelists
    config%n_sets = 0
    do
        ! Reset to invalid values
        t_1_uv = -9999.0_dpreal; t_1_rv = -9999.0_dpreal; t_1_av = -9999.0_dpreal
        t_1_uw = -9999.0_dpreal; t_1_rw = -9999.0_dpreal; t_1_aw = -9999.0_dpreal
        t_1_wd = -9999.0_dpreal; t_1_rwd = -9999.0_dpreal; t_1_awd = -9999.0_dpreal
        x_1_uv = -9999.0_dpreal; x_1_rv = -9999.0_dpreal; x_1_av = -9999.0_dpreal
        x_1_uw = -9999.0_dpreal; x_1_rw = -9999.0_dpreal; x_1_aw = -9999.0_dpreal
        x_1_wd = -9999.0_dpreal; x_1_rwd = -9999.0_dpreal; x_1_awd = -9999.0_dpreal
        b_1_uv = -9999.0_dpreal; b_1_wd = -9999.0_dpreal; p_1_uv = -9999.0_dpreal

        read(iunit, nml=PARAM_SET, iostat=ios)
        if (ios /= 0) exit

        ! Check if any parameter was set
        has_value = .false.
        if (t_1_uv > -9998.0_dpreal) has_value = .true.
        if (t_1_wd > -9998.0_dpreal) has_value = .true.
        if (x_1_uv > -9998.0_dpreal) has_value = .true.
        if (x_1_wd > -9998.0_dpreal) has_value = .true.

        if (.not. has_value) exit

        config%n_sets = config%n_sets + 1
        is = config%n_sets

        if (is > MAX_PARAM_SETS) then
            ierr = -3
            write(*,'(A,I5)') 'ERROR: Too many parameter sets, max=', MAX_PARAM_SETS
            exit
        end if

        ! Store parameters that were set
        config%sets(is)%n_params = 0

        ! Proton potential parameters
        if (t_1_uv > -9998.0_dpreal) call add_param(config%sets(is), 't_1_uv', t_1_uv)
        if (t_1_rv > -9998.0_dpreal) call add_param(config%sets(is), 't_1_rv', t_1_rv)
        if (t_1_av > -9998.0_dpreal) call add_param(config%sets(is), 't_1_av', t_1_av)
        if (t_1_uw > -9998.0_dpreal) call add_param(config%sets(is), 't_1_uw', t_1_uw)
        if (t_1_rw > -9998.0_dpreal) call add_param(config%sets(is), 't_1_rw', t_1_rw)
        if (t_1_aw > -9998.0_dpreal) call add_param(config%sets(is), 't_1_aw', t_1_aw)
        if (t_1_wd > -9998.0_dpreal) call add_param(config%sets(is), 't_1_wd', t_1_wd)
        if (t_1_rwd > -9998.0_dpreal) call add_param(config%sets(is), 't_1_rwd', t_1_rwd)
        if (t_1_awd > -9998.0_dpreal) call add_param(config%sets(is), 't_1_awd', t_1_awd)
        ! Neutron potential parameters
        if (x_1_uv > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_uv', x_1_uv)
        if (x_1_rv > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_rv', x_1_rv)
        if (x_1_av > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_av', x_1_av)
        if (x_1_uw > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_uw', x_1_uw)
        if (x_1_rw > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_rw', x_1_rw)
        if (x_1_aw > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_aw', x_1_aw)
        if (x_1_wd > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_wd', x_1_wd)
        if (x_1_rwd > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_rwd', x_1_rwd)
        if (x_1_awd > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_awd', x_1_awd)
        ! Other potentials
        if (b_1_uv > -9998.0_dpreal) call add_param(config%sets(is), 'b_1_uv', b_1_uv)
        if (b_1_wd > -9998.0_dpreal) call add_param(config%sets(is), 'b_1_wd', b_1_wd)
        if (p_1_uv > -9998.0_dpreal) call add_param(config%sets(is), 'p_1_uv', p_1_uv)
    end do

    close(iunit)

    ! Print configuration
    call print_predict_config(config)

end subroutine

!==============================================================================
! Read prediction configuration from stdin
!==============================================================================
subroutine read_predict_config_stdin(config, ierr)
    type(predict_config_t), intent(out) :: config
    integer, intent(out) :: ierr

    ! Namelist variables
    character(len=256) :: emulator_file, base_input
    character(len=16) :: predict_method
    character(len=32) :: output_format

    ! PARAM_SET variables - all 18 optical potential parameters
    real(dpreal) :: t_1_uv, t_1_rv, t_1_av  ! proton volume real
    real(dpreal) :: t_1_uw, t_1_rw, t_1_aw  ! proton volume imag
    real(dpreal) :: t_1_wd, t_1_rwd, t_1_awd ! proton surface imag
    real(dpreal) :: x_1_uv, x_1_rv, x_1_av  ! neutron volume real
    real(dpreal) :: x_1_uw, x_1_rw, x_1_aw  ! neutron volume imag
    real(dpreal) :: x_1_wd, x_1_rwd, x_1_awd ! neutron surface imag
    real(dpreal) :: b_1_uv, b_1_wd, p_1_uv  ! other potentials

    integer :: ios, is
    logical :: has_value

    namelist /EMULATOR_PREDICT/ emulator_file, base_input, predict_method, output_format
    namelist /PARAM_SET/ t_1_uv, t_1_rv, t_1_av, t_1_uw, t_1_rw, t_1_aw, t_1_wd, t_1_rwd, t_1_awd, &
                         x_1_uv, x_1_rv, x_1_av, x_1_uw, x_1_rw, x_1_aw, x_1_wd, x_1_rwd, x_1_awd, &
                         b_1_uv, b_1_wd, p_1_uv

    ierr = 0

    ! Initialize defaults
    emulator_file = ''
    base_input = ''
    predict_method = 'matrix'
    output_format = 'cross_section'

    ! Read EMULATOR_PREDICT namelist from stdin
    read(*, nml=EMULATOR_PREDICT, iostat=ios)
    if (ios /= 0) then
        ierr = -2
        write(*,'(A)') 'ERROR: Failed to read EMULATOR_PREDICT namelist from stdin'
        return
    end if

    config%emulator_file = emulator_file
    config%base_input = base_input
    config%predict_method = predict_method
    config%output_format = output_format

    ! Read PARAM_SET namelists
    config%n_sets = 0
    do
        ! Reset to invalid values
        t_1_uv = -9999.0_dpreal; t_1_rv = -9999.0_dpreal; t_1_av = -9999.0_dpreal
        t_1_uw = -9999.0_dpreal; t_1_rw = -9999.0_dpreal; t_1_aw = -9999.0_dpreal
        t_1_wd = -9999.0_dpreal; t_1_rwd = -9999.0_dpreal; t_1_awd = -9999.0_dpreal
        x_1_uv = -9999.0_dpreal; x_1_rv = -9999.0_dpreal; x_1_av = -9999.0_dpreal
        x_1_uw = -9999.0_dpreal; x_1_rw = -9999.0_dpreal; x_1_aw = -9999.0_dpreal
        x_1_wd = -9999.0_dpreal; x_1_rwd = -9999.0_dpreal; x_1_awd = -9999.0_dpreal
        b_1_uv = -9999.0_dpreal; b_1_wd = -9999.0_dpreal; p_1_uv = -9999.0_dpreal

        read(*, nml=PARAM_SET, iostat=ios)
        if (ios /= 0) exit

        ! Check if any parameter was set
        has_value = .false.
        if (t_1_uv > -9998.0_dpreal) has_value = .true.
        if (t_1_wd > -9998.0_dpreal) has_value = .true.
        if (x_1_uv > -9998.0_dpreal) has_value = .true.
        if (x_1_wd > -9998.0_dpreal) has_value = .true.

        if (.not. has_value) exit

        config%n_sets = config%n_sets + 1
        is = config%n_sets

        if (is > MAX_PARAM_SETS) then
            ierr = -3
            write(*,'(A,I5)') 'ERROR: Too many parameter sets, max=', MAX_PARAM_SETS
            exit
        end if

        ! Store parameters that were set
        config%sets(is)%n_params = 0

        ! Proton potential parameters
        if (t_1_uv > -9998.0_dpreal) call add_param(config%sets(is), 't_1_uv', t_1_uv)
        if (t_1_rv > -9998.0_dpreal) call add_param(config%sets(is), 't_1_rv', t_1_rv)
        if (t_1_av > -9998.0_dpreal) call add_param(config%sets(is), 't_1_av', t_1_av)
        if (t_1_uw > -9998.0_dpreal) call add_param(config%sets(is), 't_1_uw', t_1_uw)
        if (t_1_rw > -9998.0_dpreal) call add_param(config%sets(is), 't_1_rw', t_1_rw)
        if (t_1_aw > -9998.0_dpreal) call add_param(config%sets(is), 't_1_aw', t_1_aw)
        if (t_1_wd > -9998.0_dpreal) call add_param(config%sets(is), 't_1_wd', t_1_wd)
        if (t_1_rwd > -9998.0_dpreal) call add_param(config%sets(is), 't_1_rwd', t_1_rwd)
        if (t_1_awd > -9998.0_dpreal) call add_param(config%sets(is), 't_1_awd', t_1_awd)
        ! Neutron potential parameters
        if (x_1_uv > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_uv', x_1_uv)
        if (x_1_rv > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_rv', x_1_rv)
        if (x_1_av > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_av', x_1_av)
        if (x_1_uw > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_uw', x_1_uw)
        if (x_1_rw > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_rw', x_1_rw)
        if (x_1_aw > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_aw', x_1_aw)
        if (x_1_wd > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_wd', x_1_wd)
        if (x_1_rwd > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_rwd', x_1_rwd)
        if (x_1_awd > -9998.0_dpreal) call add_param(config%sets(is), 'x_1_awd', x_1_awd)
        ! Other potentials
        if (b_1_uv > -9998.0_dpreal) call add_param(config%sets(is), 'b_1_uv', b_1_uv)
        if (b_1_wd > -9998.0_dpreal) call add_param(config%sets(is), 'b_1_wd', b_1_wd)
        if (p_1_uv > -9998.0_dpreal) call add_param(config%sets(is), 'p_1_uv', p_1_uv)
    end do

    ! Print configuration
    call print_predict_config(config)

end subroutine

!==============================================================================
! Helper: Add parameter to param_set
!==============================================================================
subroutine add_param(pset, name, value)
    type(param_set_t), intent(inout) :: pset
    character(len=*), intent(in) :: name
    real(dpreal), intent(in) :: value

    pset%n_params = pset%n_params + 1
    pset%names(pset%n_params) = name
    pset%values(pset%n_params) = value
end subroutine

!==============================================================================
! Print training configuration
!==============================================================================
subroutine print_train_config(config)
    type(train_config_t), intent(in) :: config
    integer :: i

    write(*,'(A)') '=============================================='
    write(*,'(A)') 'Emulator Training Configuration'
    write(*,'(A)') '=============================================='
    write(*,'(A,A)')    '  Base input:    ', trim(config%base_input)
    write(*,'(A,I6)')   '  N samples:     ', config%n_samples
    write(*,'(A,I6)')   '  N basis:       ', config%n_basis
    write(*,'(A,A)')    '  Method:        ', trim(config%method)
    write(*,'(A,ES10.2)')'  Tolerance:     ', config%tol
    write(*,'(A,A)')    '  Output file:   ', trim(config%output_file)
    write(*,'(A,L6)')   '  Multi-J mode:  ', config%multi_J
    if (config%multi_J) then
        write(*,'(A,L6)')   '  Random J:      ', config%random_J
        if (config%random_J) then
            write(*,'(A,I6)')   '  J samples/set: ', config%n_J_samples
        end if
    end if
    write(*,'(A)') ''
    write(*,'(A,I3)')   '  Number of parameters: ', config%n_params
    write(*,'(A)') '  Parameter ranges:'
    do i = 1, config%n_params
        write(*,'(A,A,A,F10.3,A,F10.3)') '    ', trim(config%params(i)%name), &
            ': [', config%params(i)%min_val, ', ', config%params(i)%max_val, ']'
    end do
    write(*,'(A)') '=============================================='

end subroutine

!==============================================================================
! Print prediction configuration
!==============================================================================
subroutine print_predict_config(config)
    type(predict_config_t), intent(in) :: config
    integer :: i, j

    write(*,'(A)') '=============================================='
    write(*,'(A)') 'Emulator Prediction Configuration'
    write(*,'(A)') '=============================================='
    write(*,'(A,A)')    '  Emulator file: ', trim(config%emulator_file)
    write(*,'(A,A)')    '  Base input:    ', trim(config%base_input)
    write(*,'(A,A)')    '  Method:        ', trim(config%predict_method)
    write(*,'(A,A)')    '  Output:        ', trim(config%output_format)
    write(*,'(A)') ''
    write(*,'(A,I5)')   '  Number of parameter sets: ', config%n_sets
    do i = 1, min(5, config%n_sets)
        write(*,'(A,I3,A)', advance='no') '    Set ', i, ': '
        do j = 1, config%sets(i)%n_params
            write(*,'(A,A,F8.3)', advance='no') trim(config%sets(i)%names(j)), '=', &
                config%sets(i)%values(j)
            if (j < config%sets(i)%n_params) write(*,'(A)', advance='no') ', '
        end do
        write(*,*)
    end do
    if (config%n_sets > 5) then
        write(*,'(A,I5,A)') '    ... (', config%n_sets - 5, ' more sets)'
    end if
    write(*,'(A)') '=============================================='

end subroutine

end module emulator_io
