!==============================================================================
! Test Multi-J Emulator Prediction Accuracy
!
! Uses matrix reconstruction to compare predictions vs exact DBMM
!==============================================================================
program test_multiJ_prediction

use precision
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
use param_mapping

implicit none

type(emulator_t) :: emu
type(nch_tilde) :: atilde
type(nch_tilde), allocatable :: atilde_J(:)

integer :: ierr, iJ, ich, j, itest, i
character(len=256) :: emulator_file, base_input

! Test parameters
integer, parameter :: n_test = 3
real(dpreal) :: test_params(2, n_test)
real(dpreal) :: param_values(2)

! DBMM variables
real(dpreal) :: mu_val, ecm_val, jtot
real(dpreal), allocatable :: k_ch(:), eta_val(:), E_thresh(:)
real(dpreal), allocatable :: k_ch_all(:,:), eta_all(:,:)
integer, allocatable :: l_ch_arr(:), ch_open_arr(:), l_ch_all(:,:)
complex(dpreal), allocatable :: u_couple_arr(:,:,:)
complex(dpreal), allocatable :: S_column(:), c_exact(:), c_pred(:), c_red(:)
complex(dpreal), allocatable :: b_full(:)

! Error tracking
real(dpreal) :: error_rel, norm_exact
real(dpreal) :: avg_error, max_error
integer :: nJ, nch_max_all, parity, ntot_J

! Timing
real(dpreal) :: t_exact_start, t_exact_end, t_pred_start, t_pred_end
real(dpreal) :: time_exact, time_pred

!==============================================================================
! Configuration
!==============================================================================
emulator_file = 'emulator_multiJ_test.dat'
base_input = '../test/test_d58Ni_cdcc.in'

! Test parameters (within training range 48-58, 6-10)
test_params(1, 1) = 53.3_dpreal   ! baseline
test_params(2, 1) = 7.8_dpreal

test_params(1, 2) = 50.0_dpreal   ! lower end
test_params(2, 2) = 8.5_dpreal

test_params(1, 3) = 55.0_dpreal   ! upper end
test_params(2, 3) = 6.5_dpreal

!==============================================================================
write(*,'(A)') '=============================================================='
write(*,'(A)') '     Multi-J Emulator Prediction Accuracy Test'
write(*,'(A)') '=============================================================='
write(*,'(A)') ''

!==============================================================================
! Step 1: Initialize physics system
!==============================================================================
write(*,'(A)') 'Step 1: Initializing physics system...'

open(unit=5, file=trim(base_input), status='old', iostat=ierr)
if (ierr /= 0) then
    write(*,'(A)') 'ERROR: Cannot open base input'
    stop 1
end if

call initialize()
if (method == 1) then
    call setup_bins()
else
    call tho_wavefunction()
end if
close(5)

call alpha_3b_new()
call phibx_bound()
call phibx_bin()

mu_val = amu * (masst * massp) / (massp + masst)
ecm_val = elab * masst / (massp + masst) + epsilon_bin(1,1)

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
write(*,'(A,I4,A,I4,A)') '  Loaded multi-J emulator with ', emu%n_samples, &
    ' training samples, ', nJ, ' J values'

!==============================================================================
! Step 3: Setup channels for all J
!==============================================================================
write(*,'(A)') ''
write(*,'(A)') 'Step 3: Setting up channels...'

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
    if (allocated(atilde_J(iJ)%ch_open)) deallocate(atilde_J(iJ)%ch_open)
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
allocate(b_full(nlag * nch_max_all))

! Allocate arrays for emulator_setup_multiJ
allocate(k_ch_all(nch_max_all, nJ))
allocate(eta_all(nch_max_all, nJ))
allocate(l_ch_all(nch_max_all, nJ))

! Precompute channel info for all J
do iJ = 1, nJ
    do ich = 1, atilde_J(iJ)%nchmax
        l_ch_all(ich, iJ) = a3b%lam(atilde_J(iJ)%alpha3b(ich))
        E_thresh(ich) = epsilon_bin(atilde_J(iJ)%n(ich), &
                                   a3b%alpha2b(atilde_J(iJ)%alpha3b(ich)))

        if ((ecm_val - E_thresh(ich)) > 0.00001_dpreal) then
            k_ch_all(ich, iJ) = sqrt(2.0_dpreal * mu_val * (ecm_val - E_thresh(ich))) / hbarc
            eta_all(ich, iJ) = zp * zt * e2 * mu_val / hbarc / hbarc / k_ch_all(ich, iJ)
        else
            k_ch_all(ich, iJ) = sqrt(2.0_dpreal * mu_val * (E_thresh(ich) - ecm_val)) / hbarc
            eta_all(ich, iJ) = 0.0_dpreal
        end if
    end do
end do

! Setup emulator for multi-J (precompute reduced kinetic+boundary matrices)
write(*,'(A)') '  Setting up emulator reduced matrices...'
call emulator_setup_multiJ(emu, atilde_J, mu_val, ecm_val, k_ch_all, eta_all, l_ch_all)
write(*,'(A)') '  Emulator setup complete'

!==============================================================================
! Step 4: Run prediction tests
!==============================================================================
write(*,'(A)') ''
write(*,'(A)') '=============================================================='
write(*,'(A)') '  Prediction Accuracy Results (Matrix Reconstruction Method)'
write(*,'(A)') '=============================================================='
write(*,'(A)') ''

do itest = 1, n_test
    param_values = test_params(:, itest)

    write(*,'(A,I2,A,F7.2,A,F6.2)') 'Test ', itest, ': t_1_uv=', param_values(1), &
        ', t_1_wd=', param_values(2)

    ! Set potential parameters
    call set_potential_param('t_1_uv', param_values(1), ierr)
    call set_potential_param('t_1_wd', param_values(2), ierr)

    ! Recompute potentials
    call potr('t', 1, zb*zt, 0.0_dpreal)
    UbA = V
    call potr('x', 1, zx*zt, 0.0_dpreal)
    UxA = V

    avg_error = 0.0_dpreal
    max_error = 0.0_dpreal
    time_exact = 0.0_dpreal
    time_pred = 0.0_dpreal

    write(*,'(A)') ''
    write(*,'(A)') '     J      nch    Error(%)   t_exact(ms)  t_pred(ms)  Speedup'
    write(*,'(A)') '   -----   -----   --------   -----------  ----------  -------'

    ! Test selected J values
    do iJ = 1, nJ, 5  ! Test every 5th J value
        call cleanup_potential_cache()

        ntot_J = emu%ntot_J(iJ)

        ! Allocate per-J reduced coefficient array
        if (allocated(c_red)) deallocate(c_red)
        allocate(c_red(emu%n_basis_J(iJ)))

        ! Setup channel info
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
        call cpu_time(t_exact_start)
        c_exact = (0.0_dpreal, 0.0_dpreal)
        b_full = (0.0_dpreal, 0.0_dpreal)
        call dbmm_coupled_channel_cc(atilde_J(iJ)%nchmax, atilde_J(iJ)%nch_open, &
                                      ch_open_arr(1:atilde_J(iJ)%nch_open), &
                                      nlag, rmax, l_ch_arr(1:atilde_J(iJ)%nchmax), &
                                      k_ch(1:atilde_J(iJ)%nchmax), &
                                      eta_val(1:atilde_J(iJ)%nchmax), &
                                      u_couple_arr(1:nlag, 1:atilde_J(iJ)%nchmax, 1:atilde_J(iJ)%nchmax), &
                                      mu_val, ecm_val, E_thresh(1:atilde_J(iJ)%nchmax), &
                                      1, S_column(1:atilde_J(iJ)%nch_open), &
                                      b_out=b_full(1:ntot_J), &
                                      c_out=c_exact(1:ntot_J))
        call cpu_time(t_exact_end)

        ! Emulator prediction using matrix reconstruction
        call cpu_time(t_pred_start)
        c_pred = (0.0_dpreal, 0.0_dpreal)
        c_red = (0.0_dpreal, 0.0_dpreal)
        call emulator_predict_J(emu, iJ, &
                                u_couple_arr(1:nlag, 1:atilde_J(iJ)%nchmax, 1:atilde_J(iJ)%nchmax), &
                                b_full(1:ntot_J), c_red, c_pred(1:ntot_J))
        call cpu_time(t_pred_end)

        ! Compute relative error
        norm_exact = sqrt(sum(abs(c_exact(1:ntot_J))**2))
        error_rel = sqrt(sum(abs(c_pred(1:ntot_J) - c_exact(1:ntot_J))**2)) / (norm_exact + 1.0d-30)

        avg_error = avg_error + error_rel
        if (error_rel > max_error) max_error = error_rel
        time_exact = time_exact + (t_exact_end - t_exact_start)
        time_pred = time_pred + (t_pred_end - t_pred_start)

        write(*,'(F8.1,I8,F11.4,F13.2,F12.4,F9.1,A)') &
            emu%J_values(iJ), atilde_J(iJ)%nchmax, &
            error_rel * 100.0_dpreal, &
            (t_exact_end - t_exact_start) * 1000.0_dpreal, &
            (t_pred_end - t_pred_start) * 1000.0_dpreal, &
            (t_exact_end - t_exact_start) / (t_pred_end - t_pred_start + 1.0d-10), 'x'
    end do

    write(*,'(A)') ''
    write(*,'(A,F8.4,A)') '  Average error: ', avg_error / real((nJ+4)/5, dpreal) * 100.0_dpreal, '%'
    write(*,'(A,F8.4,A)') '  Max error:     ', max_error * 100.0_dpreal, '%'
    write(*,'(A,F8.1,A)') '  Total exact time:   ', time_exact * 1000.0_dpreal, ' ms'
    write(*,'(A,F8.4,A)') '  Total predict time: ', time_pred * 1000.0_dpreal, ' ms'
    write(*,'(A,F8.1,A)') '  Overall speedup:    ', time_exact / (time_pred + 1.0d-10), 'x'
    write(*,'(A)') ''
end do

write(*,'(A)') '=============================================================='
write(*,'(A)') '  Test Complete'
write(*,'(A)') '=============================================================='

! Cleanup
call emulator_free(emu)

end program test_multiJ_prediction
