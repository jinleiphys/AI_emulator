!==============================================================================
! Parameter Mapping Module
! Maps parameter names to potential array indices in pot_cdcc
!==============================================================================
module param_mapping

use precision
use pot_cdcc
implicit none

contains

!==============================================================================
! Find potential index by kp1 and kp2
! Returns 0 if not found
!==============================================================================
function find_potential_index(kp1, kp2) result(idx)
    character(len=*), intent(in) :: kp1
    integer, intent(in) :: kp2
    integer :: idx

    integer :: i

    idx = 0
    do i = 1, 10
        if (nkp1(i) == kp1(1:1) .and. nkp2(i) == kp2) then
            idx = i
            return
        end if
    end do

end function

!==============================================================================
! Set a single potential parameter by name
! name format: 'kp1_kp2_param' (e.g., 't_1_uv')
!==============================================================================
subroutine set_potential_param(name, value, ierr)
    character(len=*), intent(in) :: name
    real(dpreal), intent(in) :: value
    integer, intent(out) :: ierr

    character(len=4) :: kp1
    integer :: kp2
    character(len=8) :: param
    integer :: idx

    ierr = 0

    ! Parse parameter name
    call parse_param_name_internal(name, kp1, kp2, param, ierr)
    if (ierr /= 0) then
        write(*,'(A)') 'ERROR: Invalid parameter name: '//trim(name)
        return
    end if

    ! Find potential index
    idx = find_potential_index(kp1, kp2)
    if (idx == 0) then
        write(*,'(A)') 'ERROR: Potential not found: '//trim(kp1)//'_'
        write(*,'(A,I2)') '  kp2 = ', kp2
        ierr = -1
        return
    end if

    ! Set the parameter value
    select case (trim(param))
    case ('uv')
        nuv(idx) = value
    case ('rv')
        nrv(idx) = value
    case ('av')
        nav(idx) = value
    case ('uw')
        nuw(idx) = value
    case ('rw')
        nrw(idx) = value
    case ('aw')
        naw(idx) = value
    case ('wd')
        nwd(idx) = value
    case ('rwd')
        nrwd(idx) = value
    case ('awd')
        nawd(idx) = value
    case ('vd')
        nvd(idx) = value
    case ('rvd')
        nrvd(idx) = value
    case ('avd')
        navd(idx) = value
    case ('rc')
        nrc(idx) = value
    case default
        write(*,'(A)') 'ERROR: Unknown parameter: '//trim(param)
        ierr = -2
    end select

end subroutine

!==============================================================================
! Get a single potential parameter by name
!==============================================================================
function get_potential_param(name, ierr) result(value)
    character(len=*), intent(in) :: name
    integer, intent(out) :: ierr
    real(dpreal) :: value

    character(len=4) :: kp1
    integer :: kp2
    character(len=8) :: param
    integer :: idx

    ierr = 0
    value = 0.0_dpreal

    ! Parse parameter name
    call parse_param_name_internal(name, kp1, kp2, param, ierr)
    if (ierr /= 0) return

    ! Find potential index
    idx = find_potential_index(kp1, kp2)
    if (idx == 0) then
        ierr = -1
        return
    end if

    ! Get the parameter value
    select case (trim(param))
    case ('uv')
        value = nuv(idx)
    case ('rv')
        value = nrv(idx)
    case ('av')
        value = nav(idx)
    case ('uw')
        value = nuw(idx)
    case ('rw')
        value = nrw(idx)
    case ('aw')
        value = naw(idx)
    case ('wd')
        value = nwd(idx)
    case ('rwd')
        value = nrwd(idx)
    case ('awd')
        value = nawd(idx)
    case ('vd')
        value = nvd(idx)
    case ('rvd')
        value = nrvd(idx)
    case ('avd')
        value = navd(idx)
    case ('rc')
        value = nrc(idx)
    case default
        ierr = -2
    end select

end function

!==============================================================================
! Set multiple parameters from arrays
!==============================================================================
subroutine set_potential_params(n_params, names, values, ierr)
    integer, intent(in) :: n_params
    character(len=32), intent(in) :: names(n_params)
    real(dpreal), intent(in) :: values(n_params)
    integer, intent(out) :: ierr

    integer :: i, err

    ierr = 0
    do i = 1, n_params
        call set_potential_param(names(i), values(i), err)
        if (err /= 0) then
            ierr = err
            write(*,'(A,I3,A,A)') 'ERROR at parameter ', i, ': ', trim(names(i))
        end if
    end do

end subroutine

!==============================================================================
! Save current potential parameters to arrays
!==============================================================================
subroutine save_potential_params(n_params, names, values)
    integer, intent(in) :: n_params
    character(len=32), intent(in) :: names(n_params)
    real(dpreal), intent(out) :: values(n_params)

    integer :: i, err

    do i = 1, n_params
        values(i) = get_potential_param(names(i), err)
    end do

end subroutine

!==============================================================================
! Internal: Parse parameter name
!==============================================================================
subroutine parse_param_name_internal(name, kp1, kp2, param, ierr)
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
! Print all potential parameters for debugging
!==============================================================================
subroutine print_potential_params()
    integer :: i

    write(*,'(A)') '=============================================='
    write(*,'(A)') 'Current Potential Parameters'
    write(*,'(A)') '=============================================='
    write(*,'(A)') ' idx  kp1  kp2    uv       rv       av       wd       rwd      awd'
    write(*,'(A)') '----------------------------------------------------------------------'

    do i = 1, 10
        if (nkp1(i) /= ' ' .and. nkp1(i) /= '') then
            write(*,'(I3,3X,A1,3X,I2,6F9.3)') i, nkp1(i), nkp2(i), &
                nuv(i), nrv(i), nav(i), nwd(i), nrwd(i), nawd(i)
        end if
    end do
    write(*,'(A)') '=============================================='

end subroutine

end module param_mapping
