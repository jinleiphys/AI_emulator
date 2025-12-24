!==============================================================================
! Test Program for Multi-J Emulator Accuracy
!
! Compares emulator predictions with exact DBMM calculations
!==============================================================================
program test_multiJ_accuracy

use precision
use emulator_io
use param_mapping
use dbmm_emulator
use dbmm

! Cookie modules
use input
use channels_bin
use bin
use tho
use systems
use pot_cdcc
use mesh
use constants
use schennu
use ucouple
use ucouple_cached, only: cleanup_potential_cache

implicit none

type(emulator_t) :: emu
integer :: ierr, i, j, iJ, ich, itest
character(len=256) :: emulator_file, base_input

! Test parameters
integer, parameter :: n_test = 5
real(dpreal) :: test_params(2, n_test)
real(dpreal) :: param_values(2)

! CDCC variables
type(nch_tilde) :: atilde
type(nch_tilde), allocatable :: atilde_J(:)
real(dpreal) :: mu_val, ecm_val, jtot
real(dpreal), allocatable :: k_ch(:), eta_val(:), E_thresh(:)
integer, allocatable :: l_ch_arr(:), ch_open_arr(:)
complex(dpreal), allocatable :: u_couple_arr(:,:,:)
complex(dpreal), allocatable :: S_column(:), c_exact(:), c_pred(:)

! Error tracking
real(dpreal) :: error_J, max_error_J, avg_error_J
real(dpreal) :: error_total, max_error_total
real(dpreal), allocatable :: errors_per_J(:)
integer, allocatable :: samples_per_J(:)
integer :: nJ, nch_max_all, parity

! Timing
real(dpreal) :: t_start, t_end, t_exact, t_pred

!==============================================================================
! Configuration
!==============================================================================
emulator_file = 'emulator_randomJ_test.dat'
base_input = '../test/test_d58Ni_cdcc.in'

! Test parameters (within training range)
test_params(1, 1) = 50.0_dpreal   ! t_1_uv
test_params(2, 1) = 7.0_dpreal    ! t_1_wd

test_params(1, 2) = 53.3_dpreal   ! baseline
test_params(2, 2) = 7.8_dpreal

test_params(1, 3) = 55.0_dpreal
test_params(2, 3) = 8.5_dpreal

test_params(1, 4) = 48.5_dpreal
test_params(2, 4) = 9.0_dpreal

test_params(1, 5) = 57.0_dpreal
test_params(2, 5) = 6.5_dpreal

!==============================================================================
write(*,'(A)') '=============================================================='
write(*,'(A)') '     Multi-J Emulator Accuracy Test'
write(*,'(A)') '=============================================================='
write(*,'(A)') ''

!==============================================================================
! Step 1: Initialize physics system
!==============================================================================
write(*,'(A)') 'Step 1: Initializing physics system...'
call read_input(trim(base_input))
call systems_init()
call bin_init()

mu_val = m_p * m_t / (m_p + m_t)
ecm_val = Elab * m_t / (m_p + m_t)

write(*,'(A,F8.3,A)') '  Ecm = ', ecm_val, ' MeV'

!==============================================================================
! Step 2: Load emulator
!==============================================================================
write(*,'(A)') ''
write(*,'(A)') 'Step 2: Loading emulator...'
call emulator_load(emu, trim(emulator_file))

if (.not. emu%multi_J) then
    write(*,'(A)') 'ERROR: Emulator is not multi-J mode!'
    stop 1
end if

nJ = emu%nJ
write(*,'(A,I4,A)') '  Loaded multi-J emulator with ', nJ, ' J values'

!==============================================================================
! Step 3: Setup channel information for all J
!==============================================================================
write(*,'(A)') ''
write(*,'(A)') 'Step 3: Setting up channels for all J values...'

allocate(atilde_J(nJ))
nch_max_all = 0
parity = 1

do iJ = 1, nJ
    jtot = emu%J_values(iJ)
    call alphatildeindex(jtot, parity, atilde_J(iJ))

    ! Compute nch_open
    atilde_J(iJ)%nch_open = 0
    do ich = 1, atilde_J(iJ)%nchmax
        if ((ecm_val - epsilon_bin(atilde_J(iJ)%n(ich), &
             a3b%alpha2b(atilde_J(iJ)%alpha3b(ich)))) > -0.00001_dpreal) then
            atilde_J(iJ)%nch_open = atilde_J(iJ)%nch_open + 1
        end if
    end do

    ! Allocate ch_open array
    if (atilde_J(iJ)%nch_open > 0) then
        allocate(atilde_J(iJ)%ch_open(atilde_J(iJ)%nch_open))
        j = 0
        do ich = 1, atilde_J(iJ)%nchmax
            if ((ecm_val - epsilon_bin(atilde_J(iJ)%n(ich), &
                 a3b%alpha2b(atilde_J(iJ)%alpha3b(ich)))) > -0.00001_dpreal) then
                j = j + 1
                atilde_J(iJ)%ch_open(j) = ich
            end if
        end do
    end if

    if (atilde_J(iJ)%nchmax > nch_max_all) nch_max_all = atilde_J(iJ)%nchmax
end do

write(*,'(A,I4)') '  Max channels: ', nch_max_all

! Allocate working arrays
allocate(k_ch(nch_max_all))
allocate(eta_val(nch_max_all))
allocate(E_thresh(nch_max_all))
allocate(l_ch_arr(nch_max_all))
allocate(ch_open_arr(nch_max_all))
allocate(u_couple_arr(nlag, nch_max_all, nch_max_all))
allocate(S_column(nch_max_all))
allocate(c_exact(nlag * nch_max_all))
allocate(c_pred(nlag * nch_max_all))
allocate(errors_per_J(nJ))
allocate(samples_per_J(nJ))

errors_per_J = 0.0_dpreal
samples_per_J = 0

!==============================================================================
! Step 4: Run accuracy tests
!==============================================================================
write(*,'(A)') ''
write(*,'(A)') 'Step 4: Running accuracy tests...'
write(*,'(A)') ''
write(*,'(A)') '  Testing ', n_test, ' parameter sets across ', nJ, ' J values'
write(*,'(A)') ''

max_error_total = 0.0_dpreal

do itest = 1, n_test
    param_values = test_params(:, itest)

    write(*,'(A,I2,A,F7.2,A,F6.2,A)') '  Test ', itest, ': t_1_uv=', param_values(1), &
        ', t_1_wd=', param_values(2), ''

    ! Set potential parameters
    call set_potential_param('t_1_uv', param_values(1), ierr)
    call set_potential_param('t_1_wd', param_values(2), ierr)

    ! Recompute potentials
    call potr('t', 1, zb*zt, 0.0_dpreal)
    UbA = V
    call potr('x', 1, zx*zt, 0.0_dpreal)
    UxA = V

    max_error_J = 0.0_dpreal
    avg_error_J = 0.0_dpreal

    ! Test each J value
    do iJ = 1, nJ
        ! Skip J values with no basis (no training data)
        if (emu%n_basis_J(iJ) == 0) cycle

        call cleanup_potential_cache()

        ! Setup channel info for this J
        do ich = 1, atilde_J(iJ)%nchmax
            l_ch_arr(ich) = a3b%lam(atilde_J(iJ)%alpha3b(ich))
            E_thresh(ich) = epsilon_bin(atilde_J(iJ)%n(ich), &
                                       a3b%alpha2b(atilde_J(iJ)%alpha3b(ich)))

            if ((ecm_val - E_thresh(ich)) > 0.00001_dpreal) then
                k_ch(ich) = sqrt(2.0_dpreal * mu_val * (ecm_val - E_thresh(ich))) / hbarc
                eta_val(ich) = zp * zt * e2 * mu_val / hbarc / hbarc / k_ch(ich)
            else
                k_ch(ich) = sqrt(2.0_dpreal * mu_val * (E_thresh(ich) - ecm_val)) / hbarc
                eta_val(ich) = 0.0_dpreal
            end if
        end do

        do ich = 1, atilde_J(iJ)%nch_open
            ch_open_arr(ich) = ich
        end do

        ! Compute coupling potential
        u_couple_arr = (0.0_dpreal, 0.0_dpreal)
        call u_alpha_tilde(u_couple_arr(1:nlag, 1:atilde_J(iJ)%nchmax, 1:atilde_J(iJ)%nchmax), &
                          atilde_J(iJ))

        ! Exact DBMM calculation
        call cpu_time(t_start)
        c_exact = (0.0_dpreal, 0.0_dpreal)
        call dbmm_coupled_channel_cc(atilde_J(iJ)%nchmax, atilde_J(iJ)%nch_open, &
                                      ch_open_arr(1:atilde_J(iJ)%nch_open), &
                                      nlag, rmax, l_ch_arr(1:atilde_J(iJ)%nchmax), &
                                      k_ch(1:atilde_J(iJ)%nchmax), &
                                      eta_val(1:atilde_J(iJ)%nchmax), &
                                      u_couple_arr(1:nlag, 1:atilde_J(iJ)%nchmax, 1:atilde_J(iJ)%nchmax), &
                                      mu_val, ecm_val, E_thresh(1:atilde_J(iJ)%nchmax), &
                                      1, S_column(1:atilde_J(iJ)%nch_open), &
                                      c_out=c_exact(1:nlag*atilde_J(iJ)%nchmax))
        call cpu_time(t_end)
        t_exact = t_end - t_start

        ! Emulator prediction
        call cpu_time(t_start)
        c_pred = (0.0_dpreal, 0.0_dpreal)
        call emulator_predict_J(emu, iJ, param_values, c_pred(1:nlag*atilde_J(iJ)%nchmax))
        call cpu_time(t_end)
        t_pred = t_end - t_start

        ! Compute relative error
        error_J = sqrt(sum(abs(c_pred(1:nlag*atilde_J(iJ)%nchmax) - &
                              c_exact(1:nlag*atilde_J(iJ)%nchmax))**2)) / &
                  sqrt(sum(abs(c_exact(1:nlag*atilde_J(iJ)%nchmax))**2) + 1.0d-30)

        errors_per_J(iJ) = errors_per_J(iJ) + error_J
        samples_per_J(iJ) = samples_per_J(iJ) + 1

        if (error_J > max_error_J) max_error_J = error_J
        avg_error_J = avg_error_J + error_J

        ! Print details for first few J values
        if (iJ <= 3 .or. iJ == nJ .or. error_J > 0.1_dpreal) then
            write(*,'(A,F5.1,A,I3,A,ES10.2,A,F8.4,A,F8.4,A)') &
                '    J=', emu%J_values(iJ), ' (basis=', emu%n_basis_J(iJ), &
                '): error=', error_J, ', t_exact=', t_exact*1000, 'ms, t_pred=', t_pred*1000, 'ms'
        end if
    end do

    avg_error_J = avg_error_J / real(nJ, dpreal)
    if (max_error_J > max_error_total) max_error_total = max_error_J

    write(*,'(A,ES10.2,A,ES10.2)') '    => Avg error: ', avg_error_J, ', Max error: ', max_error_J
    write(*,'(A)') ''
end do

!==============================================================================
! Step 5: Summary statistics
!==============================================================================
write(*,'(A)') '=============================================================='
write(*,'(A)') '                    ACCURACY SUMMARY'
write(*,'(A)') '=============================================================='
write(*,'(A)') ''
write(*,'(A)') '  Per-J average errors:'
write(*,'(A)') '  ----------------------'

do iJ = 1, nJ
    if (samples_per_J(iJ) > 0) then
        errors_per_J(iJ) = errors_per_J(iJ) / real(samples_per_J(iJ), dpreal)
        write(*,'(A,F5.1,A,I3,A,ES10.2)') '    J=', emu%J_values(iJ), &
            ' (basis=', emu%n_basis_J(iJ), '): ', errors_per_J(iJ)
    end if
end do

write(*,'(A)') ''
write(*,'(A,ES10.2)') '  Maximum error across all tests: ', max_error_total
write(*,'(A)') ''

! Identify problematic J values
write(*,'(A)') '  J values with high error (>1%):'
do iJ = 1, nJ
    if (errors_per_J(iJ) > 0.01_dpreal) then
        write(*,'(A,F5.1,A,I3,A,ES10.2)') '    J=', emu%J_values(iJ), &
            ' (basis=', emu%n_basis_J(iJ), '): ', errors_per_J(iJ)
    end if
end do

write(*,'(A)') ''
write(*,'(A)') '=============================================================='

! Cleanup
call emulator_free(emu)

end program test_multiJ_accuracy
