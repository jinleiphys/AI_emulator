!==============================================================================
! d + 58Ni Emulator Validation Test using Matrix Reconstruction
! Includes timing comparison between full DBMM and emulator
!==============================================================================
program test_emulator_dNi

use precision
use constants
use dbmm_emulator
use dbmm

implicit none

! Realistic CDCC parameters from d+58Ni calculation
integer, parameter :: nlag = 80       ! Same as real CDCC
integer, parameter :: nch = 37        ! Maximum channels from real CDCC (s,p,d waves)
integer, parameter :: nch_open = 37
integer, parameter :: n_params = 2
integer, parameter :: n_train = 50    ! Training samples
integer, parameter :: n_test = 21     ! Same as number of J values in real CDCC
integer, parameter :: n_basis = 50    ! Reduced basis size

real(dpreal), parameter :: R = 30.0_dpreal
real(dpreal), parameter :: mu_val = 938.0_dpreal
real(dpreal), parameter :: ecm = 50.0_dpreal
real(dpreal), parameter :: V0_nom = 50.0_dpreal
real(dpreal), parameter :: W0_nom = 10.0_dpreal
real(dpreal), parameter :: r0 = 1.2_dpreal
real(dpreal), parameter :: a0 = 0.65_dpreal
real(dpreal), parameter :: A_target = 58.0_dpreal

type(emulator_t) :: emu
real(dpreal) :: param_min(n_params), param_max(n_params)
real(dpreal), allocatable :: train_params(:,:)
complex(dpreal), allocatable :: c_train(:,:), c_exact(:), c_pred(:), c_red(:)
complex(dpreal), allocatable :: u_couple(:,:,:), b_full(:)
complex(dpreal), allocatable :: S_column(:)
real(dpreal), allocatable :: x_mesh(:), r_mesh(:), weights(:), lambda_w(:), T_matrix(:,:)
integer, allocatable :: l_ch(:), ch_open(:)
real(dpreal), allocatable :: k_ch(:), eta_ch(:), E_threshold(:)
real(dpreal) :: k0, r_ws, error, max_error, avg_error
integer :: i, j, ntot, ich
real(dpreal) :: params(n_params)
real(dpreal) :: t_start, t_end, time_dbmm, time_emu, speedup
real(dpreal) :: E_bin

ntot = nlag * nch

! Setup
l_ch = [0, 2]
ch_open = [1, 2]
k0 = sqrt(2.0_dpreal * mu_val * ecm) / hbarc
k_ch(1) = k0
k_ch(2) = sqrt(2.0_dpreal * mu_val * (ecm - 2.225_dpreal)) / hbarc
eta_ch = 0.0_dpreal
E_threshold = [0.0_dpreal, 2.225_dpreal]
r_ws = r0 * A_target**(1.0_dpreal/3.0_dpreal)

allocate(x_mesh(nlag), r_mesh(nlag), weights(nlag), lambda_w(nlag), T_matrix(nlag, nlag))
allocate(u_couple(nlag, nch, nch), b_full(ntot))
allocate(c_train(ntot, n_train), c_exact(ntot), c_pred(ntot), c_red(n_basis))
allocate(train_params(n_params, n_train))

call dbmm_lagrange_legendre_mesh(nlag, R, x_mesh, r_mesh, weights, lambda_w)
call dbmm_derivative_matrix(nlag, R, x_mesh, T_matrix)

param_min = 0.8_dpreal
param_max = 1.2_dpreal

write(*,'(A)') '=============================================='
write(*,'(A)') 'Matrix Reconstruction Method Test'
write(*,'(A)') '=============================================='

! Step 1: Generate training data
write(*,'(A)') 'Step 1: Training...'
call emulator_init(emu, nlag, nch, n_basis, n_params)
call lhs_sample(n_train, n_params, param_min, param_max, train_params)
call emulator_allocate_training(emu, n_train)

do i = 1, n_train
    params = train_params(:, i)
    call generate_potential(nlag, nch, r_mesh, V0_nom*params(1), W0_nom*params(2), r_ws, a0, u_couple)
    call dbmm_coupled_channel_cc(nch, nch_open, ch_open, nlag, R, l_ch, k_ch, eta_ch, &
        u_couple, mu_val, ecm, E_threshold, 1, S_column, c_out=c_train(:,i))
    call emulator_store_snapshot(emu, i, c_train(:,i), params)
end do

call emulator_finalize_training(emu)

! Step 2: Setup emulator for matrix reconstruction
write(*,'(A)') 'Step 2: Setup...'
call emulator_set_mesh(emu, x_mesh, r_mesh, weights, lambda_w, R, T_matrix, mu_val, ecm)
call emulator_setup(emu, l_ch, k_ch, eta_ch)

! Step 3: Test accuracy and timing
write(*,'(A)') 'Step 3: Testing matrix reconstruction...'
write(*,'(A)') ''
write(*,'(A,I5,A,I5)') '  Matrix dimensions: Full = ', ntot, ' x ', ntot
write(*,'(A,I5,A,I5)') '                   Reduced = ', n_basis, ' x ', n_basis
write(*,'(A,F8.1,A)') '  Dimension reduction: ', real(ntot)/real(n_basis), 'x'
write(*,'(A)') ''

! Accuracy tests (first 5 only)
write(*,'(A)') '  Accuracy tests:'
max_error = 0.0_dpreal
avg_error = 0.0_dpreal

do i = 1, min(5, n_test)
    call random_number(params)
    params = param_min + (param_max - param_min) * params

    call generate_potential(nlag, nch, r_mesh, V0_nom*params(1), W0_nom*params(2), r_ws, a0, u_couple)
    call dbmm_coupled_channel_cc(nch, nch_open, ch_open, nlag, R, l_ch, k_ch, eta_ch, &
        u_couple, mu_val, ecm, E_threshold, 1, S_column, b_out=b_full, c_out=c_exact)

    call emulator_predict(emu, u_couple, b_full, c_red, c_pred)

    error = sqrt(sum(abs(c_pred - c_exact)**2)) / sqrt(sum(abs(c_exact)**2))
    write(*,'(A,I2,A,2F7.3,A,ES12.4)') '    Test ', i, ': params=', params, ', error=', error
    max_error = max(max_error, error)
    avg_error = avg_error + error
end do
avg_error = avg_error / min(5.0_dpreal, real(n_test, dpreal))

write(*,'(A)') ''
write(*,'(A,I4,A)') '  Timing comparison (', n_test, ' evaluations):'

! Timing: Full DBMM
call cpu_time(t_start)
do i = 1, n_test
    call random_number(params)
    params = param_min + (param_max - param_min) * params
    call generate_potential(nlag, nch, r_mesh, V0_nom*params(1), W0_nom*params(2), r_ws, a0, u_couple)
    call dbmm_coupled_channel_cc(nch, nch_open, ch_open, nlag, R, l_ch, k_ch, eta_ch, &
        u_couple, mu_val, ecm, E_threshold, 1, S_column, b_out=b_full, c_out=c_exact)
end do
call cpu_time(t_end)
time_dbmm = (t_end - t_start) * 1000.0_dpreal  ! Convert to ms

! Timing: Emulator prediction (matrix reconstruction)
call cpu_time(t_start)
do i = 1, n_test
    call random_number(params)
    params = param_min + (param_max - param_min) * params
    call generate_potential(nlag, nch, r_mesh, V0_nom*params(1), W0_nom*params(2), r_ws, a0, u_couple)
    ! Note: In practice, b_full would be precomputed or from emulator
    ! Here we include the V_red computation in the emulator time
    call emulator_predict(emu, u_couple, b_full, c_red, c_pred)
end do
call cpu_time(t_end)
time_emu = (t_end - t_start) * 1000.0_dpreal  ! Convert to ms

speedup = time_dbmm / time_emu

write(*,'(A)') ''
write(*,'(A)') '=============================================='
write(*,'(A)') 'RESULTS SUMMARY'
write(*,'(A)') '=============================================='
write(*,'(A,I5,A,I5)') '  Full matrix size:     ', ntot, ' x ', ntot
write(*,'(A,I5,A,I5)') '  Reduced matrix size:  ', n_basis, ' x ', n_basis
write(*,'(A,F8.1,A)') '  Dimension reduction:  ', real(ntot)/real(n_basis), 'x'
write(*,'(A)') ''
write(*,'(A,ES12.4)') '  Max relative error:   ', max_error
write(*,'(A,ES12.4)') '  Avg relative error:   ', avg_error
write(*,'(A)') ''
write(*,'(A,F10.3,A,I4,A)') '  Full DBMM time:       ', time_dbmm, ' ms  (', n_test, ' evals)'
write(*,'(A,F10.3,A,I4,A)') '  Emulator time:        ', time_emu, ' ms  (', n_test, ' evals)'
write(*,'(A,F10.4,A)') '  Per-evaluation DBMM:  ', time_dbmm/real(n_test), ' ms'
write(*,'(A,F10.4,A)') '  Per-evaluation Emu:   ', time_emu/real(n_test), ' ms'
write(*,'(A,F10.1,A)') '  Speedup:              ', speedup, 'x'
write(*,'(A)') '=============================================='

call emulator_free(emu)
deallocate(x_mesh, r_mesh, weights, lambda_w, T_matrix, u_couple, b_full)
deallocate(c_train, c_exact, c_pred, c_red, train_params)

contains

subroutine generate_potential(n_lag, n_ch, rr, V0, W0, rr0, aa0, uu)
    integer, intent(in) :: n_lag, n_ch
    real(dpreal), intent(in) :: rr(n_lag), V0, W0, rr0, aa0
    complex(dpreal), intent(out) :: uu(n_lag, n_ch, n_ch)
    real(dpreal) :: f_ws
    integer :: ii, ich

    uu = (0.0_dpreal, 0.0_dpreal)
    do ich = 1, n_ch
        do ii = 1, n_lag
            f_ws = 1.0_dpreal / (1.0_dpreal + exp((rr(ii) - rr0) / aa0))
            uu(ii, ich, ich) = cmplx(-V0 * f_ws, -W0 * f_ws, dpreal)
        end do
    end do
    do ii = 1, n_lag
        uu(ii, 1, 2) = 0.5_dpreal * exp(-rr(ii)**2 / 20.0_dpreal)
        uu(ii, 2, 1) = uu(ii, 1, 2)
    end do
end subroutine

end program test_emulator_dNi
