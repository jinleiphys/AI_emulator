!==============================================================================
! Emulator Prediction Program for CDCC Calculations
!
! Usage: ./emulator_predict <config_file>
! Example: ./emulator_predict predict.in
!
! This program:
! 1. Loads a trained emulator
! 2. Reads parameter sets to predict
! 3. Uses RBF interpolation or matrix reconstruction for fast prediction
! 4. Outputs predicted cross sections
!==============================================================================
program emulator_predict

use precision
use emulator_io
use param_mapping
use dbmm_emulator, emulator_predict_matrix => emulator_predict
use dbmm

! Cookie modules (needed for physics setup)
use input
use channels_bin
use bin
use tho
use systems
use pot_cdcc
use mesh
use constants
use schennu
use ucouple, only: u_alpha_tilde, u_alpha_tilde_selective
use ucouple_cached, only: cleanup_potential_cache
use gauss, only: gauleg

implicit none

type(predict_config_t) :: config
type(emulator_t) :: emu
integer :: ierr, i, j, is, ich, iJ, k

! Prediction arrays
complex(dpreal), allocatable :: c_pred(:), c_red(:), c_exact(:), c_eim(:)
complex(dpreal), allocatable :: c_red_J(:), c_full_J(:)  ! Per-J prediction arrays
complex(dpreal), allocatable :: S_pred(:), S_exact(:), S_column(:)
complex(dpreal), allocatable :: V_interp(:)  ! EIM interpolation values
real(dpreal), allocatable :: params(:)
real(dpreal) :: error_eim
real(dpreal) :: t_eim_start, t_eim_end

! CDCC variables for matrix reconstruction
type(nch_tilde) :: atilde
type(nch_tilde), allocatable :: atilde_J(:)  ! Per-J channel info for multi-J
real(dpreal) :: mu_val, ecm_val, jtot
real(dpreal), allocatable :: k_ch(:), eta_val(:), E_thresh(:)
integer, allocatable :: l_ch_arr(:), ch_open_arr(:)
complex(dpreal), allocatable :: u_couple_arr(:,:,:), b_full(:)
! Note: kinetic_matrix, l_barrier, Bc1 removed - DBMM handles internally
real(dpreal), allocatable :: x_mesh(:), r_mesh_arr(:), weights(:), lambda_arr(:)
real(dpreal), allocatable :: T_matrix(:,:)

! Multi-J variables
integer :: nch_max_all
real(dpreal) :: sigma_el_total, sigma_el_J
real(dpreal) :: error_J, norm_exact_J
real(dpreal), allocatable :: k_ch_all(:,:), eta_all(:,:)
integer, allocatable :: l_ch_all(:,:)

! S-matrix and cross section comparison
complex(dpreal), allocatable :: S_matrix(:)           ! S-matrix from predictions
complex(dpreal) :: S11_exact, S11_matrix              ! Elastic S-matrix element
real(dpreal) :: error_S_matrix                        ! S-matrix errors (per J)
real(dpreal) :: error_S_matrix_total                  ! Accumulated S-matrix errors
real(dpreal) :: sigma_exact, sigma_matrix             ! Cross sections (accumulated over J)
real(dpreal) :: dsigma_exact, dsigma_matrix           ! Per-J contribution
real(dpreal) :: k_elastic                             ! Elastic channel wavenumber
real(dpreal) :: jbx                                   ! Spin of ground state (for cross section weight)
real(dpreal), allocatable :: x_mesh_stored(:)         ! Stored mesh for S computation

! Timing for multi-J efficiency comparison (wall time using system_clock)
integer(8) :: clock_start, clock_end, clock_rate
integer(8) :: clock_exact_start, clock_exact_end, clock_matrix_start, clock_matrix_end
real(dpreal) :: time_exact_total, time_matrix_total
real(dpreal) :: time_exact_J, time_matrix_J

! File output for plotting
character(len=256) :: wf_filename, xs_filename, smat_filename
integer :: iunit_wf, iunit_xs, iunit_smat, ir
complex(dpreal), allocatable :: psi_exact_r(:), psi_matrix_r(:)
real(dpreal), allocatable :: sigma_J_exact(:), sigma_J_matrix(:)
complex(dpreal), allocatable :: S11_J_exact(:), S11_J_matrix(:)  ! S-matrix per J

! Timing and accuracy (wall time)
integer(8) :: clock_t_start, clock_t_end, clock_pred_start, clock_pred_end
real(dpreal) :: time_total, time_per_pred
real(dpreal) :: error, sigma_el
integer :: njtot, parity

! Additional timing for EIM speedup comparison (wall time)
integer(8) :: clock_pot_start, clock_pot_end, clock_eim_start, clock_eim_end
real(dpreal) :: time_full_pot, time_eim_selective
complex(dpreal), allocatable :: c_eim_fast(:)  ! EIM with selective potential
complex(dpreal), allocatable :: V_interp_sel(:)  ! Selective potential values
real(dpreal) :: error_eim_fast

!==============================================================================
! Step 0: Read prediction configuration
!==============================================================================
write(*,'(A)') '=============================================================='
write(*,'(A)') '        DBMM Emulator Prediction Program'
write(*,'(A)') '=============================================================='
write(*,'(A)') ''

! Read config directly from stdin
call read_predict_config_stdin(config, ierr)
if (ierr /= 0) then
    write(*,'(A)') 'ERROR: Failed to read prediction configuration'
    stop 1
end if

!==============================================================================
! Step 1: Load trained emulator
!==============================================================================
write(*,'(A)') ''
write(*,'(A)') 'Step 1: Loading trained emulator...'
write(*,'(A)') '  File: '//trim(config%emulator_file)

call emulator_load(emu, trim(config%emulator_file))

if (.not. emu%is_trained) then
    write(*,'(A)') 'ERROR: Emulator not trained or failed to load'
    stop 1
end if

if (emu%multi_J) then
    write(*,'(A,I6)') '  nlag:       ', emu%nlag
    write(*,'(A,I6)') '  nJ:         ', emu%nJ
    write(*,'(A,I6)') '  nch_max:    ', emu%nch_max
    write(*,'(A,I6)') '  n_params:   ', emu%n_params
    write(*,'(A,I6)') '  n_samples:  ', emu%n_samples
    write(*,'(A)') '  Multi-J emulator loaded successfully'
else
    write(*,'(A,I6)') '  nlag:       ', emu%nlag
    write(*,'(A,I6)') '  nch:        ', emu%nch
    write(*,'(A,I6)') '  ntot:       ', emu%ntot
    write(*,'(A,I6)') '  n_basis:    ', emu%n_basis
    write(*,'(A,I6)') '  n_params:   ', emu%n_params
    write(*,'(A,I6)') '  n_samples:  ', emu%n_samples
    write(*,'(A)') '  Emulator loaded successfully'
end if

!==============================================================================
! Step 2: Initialize physics system
!==============================================================================
write(*,'(A)') ''
write(*,'(A)') 'Step 2: Initializing physics system...'
write(*,'(A)') '  Base input: '//trim(config%base_input)

open(unit=5, file=trim(config%base_input), status='old', iostat=ierr)
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

! Compute basic physics parameters
ecm_val = elab * masst / (massp + masst)
mu_val = amu * (masst * massp) / (massp + masst)
ecm_val = ecm_val + epsilon_bin(1,1)
parity = 1

write(*,'(A,F8.3,A)') '  Ecm = ', ecm_val, ' MeV'

if (emu%multi_J) then
    !========================================================================
    ! Multi-J initialization
    !========================================================================
    allocate(atilde_J(emu%nJ))
    nch_max_all = 0

    ! Setup channel info for all J values
    do iJ = 1, emu%nJ
        jtot = emu%J_values(iJ)
        ! For zero spin: parity = (-1)^J (even J -> +1, odd J -> -1)
        parity = 1 - 2 * mod(nint(jtot), 2)
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

    ! Allocate working arrays with max size
    allocate(k_ch(nch_max_all))
    allocate(eta_val(nch_max_all))
    allocate(E_thresh(nch_max_all))
    allocate(l_ch_arr(nch_max_all))
    allocate(ch_open_arr(nch_max_all))
    allocate(u_couple_arr(nlag, nch_max_all, nch_max_all))
    allocate(S_column(nch_max_all))
    allocate(c_exact(nlag * nch_max_all))
    allocate(b_full(nlag * nch_max_all))

    ! Setup mesh for emulator (needed for Matrix Reconstruction)
    allocate(x_mesh(nlag), r_mesh_arr(nlag), weights(nlag), lambda_arr(nlag))
    allocate(T_matrix(nlag, nlag))
    call dbmm_lagrange_legendre_mesh(nlag, rmax, x_mesh, r_mesh_arr, weights, lambda_arr)
    call dbmm_derivative_matrix(nlag, rmax, x_mesh, T_matrix)
    call emulator_set_mesh(emu, x_mesh, r_mesh_arr, weights, lambda_arr, rmax, T_matrix, mu_val, ecm_val)

    ! Allocate arrays for emulator_setup_multiJ
    allocate(k_ch_all(nch_max_all, emu%nJ))
    allocate(eta_all(nch_max_all, emu%nJ))
    allocate(l_ch_all(nch_max_all, emu%nJ))

    ! Precompute channel info for all J (needed for Matrix Reconstruction)
    do iJ = 1, emu%nJ
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

    ! Setup emulator for Matrix Reconstruction (precompute reduced kinetic+boundary matrices)
    write(*,'(A)') '  Setting up emulator for Matrix Reconstruction...'
    call emulator_setup_multiJ(emu, atilde_J, mu_val, ecm_val, k_ch_all, eta_all, l_ch_all)
    write(*,'(A)') '  Emulator setup complete'

    write(*,'(A,I4,A)') '  Multi-J: ', emu%nJ, ' J values initialized'

else
    !========================================================================
    ! Single-J initialization
    !========================================================================
    njtot = nint(2.0_dpreal * a3b%J(1))
    jtot = njtot / 2.0_dpreal
    ! For zero spin: parity = (-1)^J (even J -> +1, odd J -> -1)
    parity = 1 - 2 * mod(njtot/2, 2)
    call alphatildeindex(jtot, parity, atilde)

    ! Compute nch_open
    atilde%nch_open = 0
    do ich = 1, atilde%nchmax
        if ((ecm_val - epsilon_bin(atilde%n(ich), a3b%alpha2b(atilde%alpha3b(ich)))) > -0.00001_dpreal) then
            atilde%nch_open = atilde%nch_open + 1
        end if
    end do

    ! Allocate and set ch_open array
    if (allocated(atilde%ch_open)) deallocate(atilde%ch_open)
    if (atilde%nch_open > 0) then
        allocate(atilde%ch_open(atilde%nch_open))
        i = 0
        do ich = 1, atilde%nchmax
            if ((ecm_val - epsilon_bin(atilde%n(ich), a3b%alpha2b(atilde%alpha3b(ich)))) > -0.00001_dpreal) then
                i = i + 1
                atilde%ch_open(i) = ich
            end if
        end do
    end if

    ! Allocate arrays
    allocate(k_ch(atilde%nchmax))
    allocate(eta_val(atilde%nchmax))
    allocate(E_thresh(atilde%nchmax))
    allocate(l_ch_arr(atilde%nchmax))
    allocate(ch_open_arr(atilde%nch_open))
    allocate(u_couple_arr(nlag, atilde%nchmax, atilde%nchmax))
    allocate(b_full(nlag * atilde%nchmax))

    ! Setup mesh for emulator
    allocate(x_mesh(nlag), r_mesh_arr(nlag), weights(nlag), lambda_arr(nlag))
    allocate(T_matrix(nlag, nlag))

    call dbmm_lagrange_legendre_mesh(nlag, rmax, x_mesh, r_mesh_arr, weights, lambda_arr)
    call dbmm_derivative_matrix(nlag, rmax, x_mesh, T_matrix)

    ! Setup channel info
    do ich = 1, atilde%nchmax
        l_ch_arr(ich) = a3b%lam(atilde%alpha3b(ich))
        E_thresh(ich) = epsilon_bin(atilde%n(ich), a3b%alpha2b(atilde%alpha3b(ich)))
        if ((ecm_val - E_thresh(ich)) > 0.00001_dpreal) then
            k_ch(ich) = sqrt(2.0_dpreal * mu_val * (ecm_val - E_thresh(ich))) / hbarc
            eta_val(ich) = zp * zt * e2 * mu_val / hbarc / hbarc / k_ch(ich)
        else
            k_ch(ich) = sqrt(2.0_dpreal * mu_val * (E_thresh(ich) - ecm_val)) / hbarc
            eta_val(ich) = 0.0_dpreal
        end if
    end do

    ! Setup emulator for matrix reconstruction
    call emulator_set_mesh(emu, x_mesh, r_mesh_arr, weights, lambda_arr, rmax, T_matrix, mu_val, ecm_val)
    call emulator_setup(emu, l_ch_arr, k_ch, eta_val)

    ! Compute b_full ONCE
    call potr('t', 1, zb*zt, 0.0_dpreal)
    UbA = V
    call potr('x', 1, zx*zt, 0.0_dpreal)
    UxA = V
    call u_alpha_tilde(u_couple_arr, atilde)
    do ich = 1, atilde%nch_open
        ch_open_arr(ich) = ich
    end do
    allocate(S_pred(atilde%nch_open), S_exact(atilde%nch_open))
    allocate(c_exact(nlag * atilde%nchmax))
    call dbmm_coupled_channel_cc(atilde%nchmax, atilde%nch_open, ch_open_arr, &
                                  nlag, rmax, l_ch_arr, k_ch, eta_val, &
                                  u_couple_arr, mu_val, ecm_val, E_thresh, &
                                  1, S_pred, b_out=b_full, c_out=c_exact)
    deallocate(S_pred, S_exact)

    write(*,'(A)') '  Source term b_full computed'

    ! Setup EIM reduced integrals if EIM is trained
    if (emu%eim_trained) then
        write(*,'(A)') ''
        write(*,'(A)') 'Step 2b: Setting up EIM reduced integrals...'
        call eim_setup_reduced_integrals(emu)
        write(*,'(A,I4,A)') '  EIM ready with ', emu%n_pot_basis, ' potential basis functions'
    end if
end if

write(*,'(A)') '  Physics system initialized'

!==============================================================================
! Step 3: Allocate prediction arrays
!==============================================================================
allocate(params(emu%n_params))

if (emu%multi_J) then
    ! For multi-J, allocate per-J arrays with max size
    allocate(c_red_J(maxval(emu%n_basis_J)))
    allocate(c_full_J(nlag * nch_max_all))
    ! S-matrix arrays for comparison
    allocate(S_matrix(nch_max_all))
    ! Store mesh for S-matrix computation
    allocate(x_mesh_stored(nlag))
    call gauleg(nlag, 0.0_dpreal, 1.0_dpreal, x_mesh_stored, weights)

    ! Allocate arrays for file output
    allocate(psi_exact_r(nlag))
    allocate(psi_matrix_r(nlag))
    allocate(sigma_J_exact(emu%nJ))
    allocate(sigma_J_matrix(emu%nJ))
    allocate(S11_J_exact(emu%nJ))
    allocate(S11_J_matrix(emu%nJ))
else
    allocate(c_pred(emu%ntot))
    allocate(c_red(emu%n_basis))

    ! Allocate EIM arrays if EIM is trained
    if (emu%eim_trained) then
        allocate(c_eim(emu%ntot))
        allocate(V_interp(emu%n_pot_interp))
        allocate(c_eim_fast(emu%ntot))
        allocate(V_interp_sel(emu%n_pot_interp))
    end if
end if

!==============================================================================
! Step 4: Run predictions
!==============================================================================
write(*,'(A)') ''
write(*,'(A)') 'Step 3: Running predictions...'
write(*,'(A)') ''

call system_clock(clock_t_start, clock_rate)

do is = 1, config%n_sets
    call system_clock(clock_pred_start)

    ! Extract parameters from param_set
    do i = 1, config%sets(is)%n_params
        params(i) = config%sets(is)%values(i)
    end do

    ! Set potential parameters
    do i = 1, config%sets(is)%n_params
        call set_potential_param(config%sets(is)%names(i), config%sets(is)%values(i), ierr)
    end do

    ! Recompute potentials
    call potr('t', 1, zb*zt, 0.0_dpreal)
    UbA = V
    call potr('x', 1, zx*zt, 0.0_dpreal)
    UxA = V

    ! Clear potential cache
    call cleanup_potential_cache()

    ! Output prediction parameters
    write(*,'(A,I4,A)', advance='no') '  Pred ', is, ': '
    do i = 1, min(config%sets(is)%n_params, 2)
        write(*,'(A,A,A,F7.2)', advance='no') trim(config%sets(is)%names(i)), '=', &
            '', config%sets(is)%values(i)
        if (i < min(config%sets(is)%n_params, 2)) write(*,'(A)', advance='no') ', '
    end do
    write(*,'(A)') ''

    if (emu%multi_J) then
        !====================================================================
        ! Multi-J prediction
        !====================================================================
        error = 0.0_dpreal
        time_exact_total = 0.0_dpreal
        time_matrix_total = 0.0_dpreal
        ! Initialize S-matrix and cross section accumulators
        error_S_matrix_total = 0.0_dpreal
        sigma_exact = 0.0_dpreal
        sigma_matrix = 0.0_dpreal

        write(*,'(A)') '      J       nch    Err_c(%)     Err_S(%)     t_exact(ms)  t_matrix(ms)  speedup'
        write(*,'(A)') '    ------   -----   ----------   ----------   -----------  ------------  -------'

        ! Loop over all J values
        do iJ = 1, emu%nJ
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

            ! Compute coupling potential for this J
            u_couple_arr = (0.0_dpreal, 0.0_dpreal)
            call u_alpha_tilde(u_couple_arr(1:nlag, 1:atilde_J(iJ)%nchmax, 1:atilde_J(iJ)%nchmax), &
                              atilde_J(iJ))

            ! Compute exact DBMM solution (with timing)
            call system_clock(clock_exact_start)
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
                                          b_out=b_full(1:emu%ntot_J(iJ)), &
                                          c_out=c_exact(1:emu%ntot_J(iJ)))
            call system_clock(clock_exact_end)
            time_exact_J = real(clock_exact_end - clock_exact_start, dpreal) / real(clock_rate, dpreal)
            time_exact_total = time_exact_total + time_exact_J

            ! Get exact S-matrix
            S11_exact = S_column(1)
            norm_exact_J = sqrt(sum(abs(c_exact(1:emu%ntot_J(iJ)))**2))

            ! Matrix reconstruction (emulator_predict_J) with timing
            call system_clock(clock_matrix_start)
            c_red_J = (0.0_dpreal, 0.0_dpreal)
            c_full_J = (0.0_dpreal, 0.0_dpreal)
            call emulator_predict_J(emu, iJ, &
                                    u_couple_arr(1:nlag, 1:atilde_J(iJ)%nchmax, 1:atilde_J(iJ)%nchmax), &
                                    b_full(1:emu%ntot_J(iJ)), &
                                    c_red_J(1:emu%n_basis_J(iJ)), &
                                    c_full_J(1:emu%ntot_J(iJ)))
            call system_clock(clock_matrix_end)
            time_matrix_J = real(clock_matrix_end - clock_matrix_start, dpreal) / real(clock_rate, dpreal)
            time_matrix_total = time_matrix_total + time_matrix_J

            ! Compute matrix reconstruction error (wavefunction coefficients)
            error = error + sqrt(sum(abs(c_full_J(1:emu%ntot_J(iJ)) - c_exact(1:emu%ntot_J(iJ)))**2)) &
                     / (norm_exact_J + 1.0d-30)

            ! Compute S-matrix from Matrix reconstruction prediction
            S_matrix = (0.0_dpreal, 0.0_dpreal)
            call compute_S_from_c(c_full_J(1:emu%ntot_J(iJ)), atilde_J(iJ)%nchmax, &
                                  atilde_J(iJ)%nch_open, ch_open_arr(1:atilde_J(iJ)%nch_open), &
                                  nlag, rmax, x_mesh_stored, &
                                  l_ch_arr(1:atilde_J(iJ)%nchmax), &
                                  k_ch(1:atilde_J(iJ)%nchmax), &
                                  eta_val(1:atilde_J(iJ)%nchmax), &
                                  1, S_matrix(1:atilde_J(iJ)%nch_open))

            ! S-matrix error for Matrix reconstruction (elastic channel S11)
            S11_matrix = S_matrix(1)
            error_S_matrix = abs(S11_matrix - S11_exact) / (abs(S11_exact) + 1.0d-30)
            error_S_matrix_total = error_S_matrix_total + error_S_matrix

            ! Accumulate elastic cross section contribution: (2J+1)/(2*jbx+1) * |1 - S11|^2
            ! σ_el = (π/k²) * Σ_J (2J+1)/(2*jbx+1) * |1 - S_J|²
            ! jbx is the angular momentum of the ground state channel (from a3b%J2b)
            k_elastic = k_ch(1)  ! Ground state channel wavenumber
            jbx = a3b%J2b(atilde_J(iJ)%alpha3b(atilde_J(iJ)%ch_open(1)))  ! Ground state J2b
            dsigma_exact = (2.0_dpreal * emu%J_values(iJ) + 1.0_dpreal) / (2.0_dpreal * jbx + 1.0_dpreal) &
                           * abs(1.0_dpreal - S11_exact)**2
            dsigma_matrix = (2.0_dpreal * emu%J_values(iJ) + 1.0_dpreal) / (2.0_dpreal * jbx + 1.0_dpreal) &
                            * abs(1.0_dpreal - S11_matrix)**2
            sigma_exact = sigma_exact + dsigma_exact
            sigma_matrix = sigma_matrix + dsigma_matrix

            ! Store per-J cross section for file output (in mb, per J)
            sigma_J_exact(iJ) = pi / k_elastic**2 * dsigma_exact * 10.0_dpreal
            sigma_J_matrix(iJ) = pi / k_elastic**2 * dsigma_matrix * 10.0_dpreal

            ! Store per-J S-matrix (elastic channel S11)
            S11_J_exact(iJ) = S11_exact
            S11_J_matrix(iJ) = S11_matrix

            ! Store wavefunction for J=0 (first J value) for elastic channel
            if (iJ == 1) then
                do ir = 1, nlag
                    psi_exact_r(ir) = c_exact(ir)
                    psi_matrix_r(ir) = c_full_J(ir)
                end do
            end if

            ! Compute error for this J
            error_J = sqrt(sum(abs(c_full_J(1:emu%ntot_J(iJ)) - c_exact(1:emu%ntot_J(iJ)))**2)) &
                     / (norm_exact_J + 1.0d-30)

            ! Output per-J result
            write(*,'(F10.1,I8,F12.4,F12.4,F13.2,F14.4,F9.1)') &
                emu%J_values(iJ), atilde_J(iJ)%nchmax, &
                error_J * 100.0_dpreal, &
                error_S_matrix * 100.0_dpreal, &
                time_exact_J * 1000.0_dpreal, &
                time_matrix_J * 1000.0_dpreal, &
                time_exact_J / (time_matrix_J + 1.0d-10)

            call cleanup_potential_cache()
        end do

        ! Average error over all J
        error = error / real(emu%nJ, dpreal)
        error_S_matrix_total = error_S_matrix_total / real(emu%nJ, dpreal)

        ! Convert cross section sum to actual cross section: σ = (π/k²) * Σ(...)
        ! Units: 10 fm² = 1000 mb, so multiply by 10 for mb
        sigma_exact = pi / k_elastic**2 * sigma_exact * 10.0_dpreal
        sigma_matrix = pi / k_elastic**2 * sigma_matrix * 10.0_dpreal

        write(*,'(A)') ''
        write(*,'(A)') '    ============== Summary =============='
        write(*,'(A)') '    Errors (avg over J):'
        write(*,'(A,F10.4,A)') '      Wavefunction: ', error * 100.0_dpreal, '%'
        write(*,'(A,F10.4,A)') '      S-matrix:     ', error_S_matrix_total * 100.0_dpreal, '%'
        write(*,'(A)') ''
        write(*,'(A)') '    Elastic cross section (nuclear):'
        write(*,'(A,F12.4,A)') '      Exact:  ', sigma_exact, ' mb'
        write(*,'(A,F12.4,A,F8.4,A)') '      Matrix: ', sigma_matrix, ' mb  (error: ', &
            abs(sigma_matrix - sigma_exact) / (sigma_exact + 1.0d-30) * 100.0_dpreal, '%)'
        write(*,'(A)') ''
        write(*,'(A)') '    Timing (all J summed):'
        write(*,'(A,F10.2,A)') '      Exact DBMM:         ', time_exact_total * 1000.0_dpreal, ' ms'
        write(*,'(A,F10.4,A,F8.1,A)') '      Matrix Reconstruct: ', time_matrix_total * 1000.0_dpreal, ' ms  (', &
            time_exact_total / (time_matrix_total + 1.0d-10), 'x speedup)'

        !====================================================================
        ! Write output files for plotting
        !====================================================================
        ! Write cross section vs J file
        write(xs_filename, '(A,I0,A)') 'cross_section_pred', is, '.dat'
        iunit_xs = 101
        open(unit=iunit_xs, file=trim(xs_filename), status='replace', action='write')
        write(iunit_xs, '(A)') '# Cross section vs J for plotting'
        write(iunit_xs, '(A)') '# Columns: J, sigma_exact(mb), sigma_matrix(mb)'
        do iJ = 1, emu%nJ
            write(iunit_xs, '(F8.1,2ES16.8)') emu%J_values(iJ), &
                sigma_J_exact(iJ), sigma_J_matrix(iJ)
        end do
        close(iunit_xs)
        write(*,'(A,A)') '    Cross section file: ', trim(xs_filename)

        ! Write wavefunction file (J=0, elastic channel)
        write(wf_filename, '(A,I0,A)') 'wavefunction_pred', is, '.dat'
        iunit_wf = 102
        open(unit=iunit_wf, file=trim(wf_filename), status='replace', action='write')
        write(iunit_wf, '(A)') '# Elastic channel wavefunction coefficients for J=0'
        write(iunit_wf, '(A)') '# Columns: r(fm), Re(c_exact), Im(c_exact), Re(c_matrix), Im(c_matrix)'
        do ir = 1, nlag
            write(iunit_wf, '(F12.6,4ES16.8)') r_mesh_arr(ir), &
                real(psi_exact_r(ir)), aimag(psi_exact_r(ir)), &
                real(psi_matrix_r(ir)), aimag(psi_matrix_r(ir))
        end do
        close(iunit_wf)
        write(*,'(A,A)') '    Wavefunction file:  ', trim(wf_filename)

        ! Write S-matrix file (elastic S11 vs J)
        write(smat_filename, '(A,I0,A)') 'smatrix_pred', is, '.dat'
        iunit_smat = 103
        open(unit=iunit_smat, file=trim(smat_filename), status='replace', action='write')
        write(iunit_smat, '(A)') '# Elastic S-matrix (S11) vs J for plotting'
        write(iunit_smat, '(A)') '# Columns: J, Re(S_exact), Im(S_exact), Re(S_matrix), Im(S_matrix)'
        do iJ = 1, emu%nJ
            write(iunit_smat, '(F8.1,4ES16.8)') emu%J_values(iJ), &
                real(S11_J_exact(iJ)), aimag(S11_J_exact(iJ)), &
                real(S11_J_matrix(iJ)), aimag(S11_J_matrix(iJ))
        end do
        close(iunit_smat)
        write(*,'(A,A)') '    S-matrix file:      ', trim(smat_filename)

    else
        !====================================================================
        ! Single-J prediction
        !====================================================================
        ! Time the full potential computation
        call system_clock(clock_pot_start)
        call u_alpha_tilde(u_couple_arr, atilde)
        call system_clock(clock_pot_end)
        time_full_pot = real(clock_pot_end - clock_pot_start, dpreal) / real(clock_rate, dpreal)

        ! Method 1: Matrix reconstruction (Galerkin projection)
        call emulator_predict_matrix(emu, u_couple_arr, b_full, c_red, c_pred)

        call system_clock(clock_pred_end)

        ! Method 2: EIM with full potential (for validation)
        if (emu%eim_trained) then
            call system_clock(clock_eim_start)
            call eim_extract_interp_values(emu, u_couple_arr, V_interp)
            call emulator_predict_eim(emu, V_interp, b_full, c_red, c_eim)
            call system_clock(clock_eim_end)
        end if

        ! Method 3: EIM with SELECTIVE potential (full speedup!)
        if (emu%eim_trained) then
            call cleanup_potential_cache()
            call system_clock(clock_eim_start)
            call u_alpha_tilde_selective(atilde, nlag, atilde%nchmax, &
                                         emu%n_pot_interp, emu%eim_indices, V_interp_sel)
            call emulator_predict_eim(emu, V_interp_sel, b_full, c_red, c_eim_fast)
            call system_clock(clock_eim_end)
            time_eim_selective = real(clock_eim_end - clock_eim_start, dpreal) / real(clock_rate, dpreal)
        end if

        ! Compute exact solution with DBMM
        allocate(S_pred(atilde%nch_open))
        call dbmm_coupled_channel_cc(atilde%nchmax, atilde%nch_open, ch_open_arr, &
                                      nlag, rmax, l_ch_arr, k_ch, eta_val, &
                                      u_couple_arr, mu_val, ecm_val, E_thresh, &
                                      1, S_pred, c_out=c_exact)
        deallocate(S_pred)

        ! Compute error for matrix reconstruction
        error = 0.0_dpreal
        do i = 1, emu%ntot
            error = error + abs(c_pred(i) - c_exact(i))**2
        end do
        error = sqrt(error / real(emu%ntot, dpreal))

        ! Compute error for EIM method (if trained)
        error_eim = 0.0_dpreal
        error_eim_fast = 0.0_dpreal
        if (emu%eim_trained) then
            do i = 1, emu%ntot
                error_eim = error_eim + abs(c_eim(i) - c_exact(i))**2
                error_eim_fast = error_eim_fast + abs(c_eim_fast(i) - c_exact(i))**2
            end do
            error_eim = sqrt(error_eim / real(emu%ntot, dpreal))
            error_eim_fast = sqrt(error_eim_fast / real(emu%ntot, dpreal))
        end if

        ! Output prediction result
        if (emu%eim_trained) then
            write(*,'(A,ES9.2,A,F7.1,A)') '    Full potential:   ', error, &
                '  (', time_full_pot*1000.0_dpreal, ' ms)'
            write(*,'(A,ES9.2,A,F7.1,A,F5.1,A)') '    EIM (selective):  ', error_eim_fast, &
                '  (', time_eim_selective*1000.0_dpreal, ' ms) => ', &
                time_full_pot/time_eim_selective, 'x speedup'
        else
            write(*,'(A,ES9.2,A,F7.1,A)') '    Matrix:', error, &
                '  (', real(clock_pred_end - clock_pred_start, dpreal) / real(clock_rate, dpreal) * 1000.0_dpreal, ' ms)'
        end if
    end if

    call system_clock(clock_pred_end)
    write(*,'(A)') ''
end do

call system_clock(clock_t_end)
time_total = real(clock_t_end - clock_t_start, dpreal) / real(clock_rate, dpreal)
if (config%n_sets > 0) then
    time_per_pred = time_total / real(config%n_sets, dpreal)
else
    time_per_pred = 0.0_dpreal
end if

!==============================================================================
! Summary
!==============================================================================
write(*,'(A)') ''
write(*,'(A)') '=============================================================='
write(*,'(A)') '                  PREDICTION COMPLETE'
write(*,'(A)') '=============================================================='
write(*,'(A,I5)')     '  Predictions made:       ', config%n_sets
write(*,'(A,F10.3,A)')'  Total prediction time:  ', time_total * 1000.0_dpreal, ' ms'
write(*,'(A,F10.3,A)')'  Time per prediction:    ', time_per_pred * 1000.0_dpreal, ' ms'
write(*,'(A)') '=============================================================='

! Cleanup
deallocate(params)
if (allocated(c_pred)) deallocate(c_pred)
if (allocated(c_red)) deallocate(c_red)
if (allocated(c_red_J)) deallocate(c_red_J)
if (allocated(c_full_J)) deallocate(c_full_J)
if (allocated(c_eim)) deallocate(c_eim)
if (allocated(V_interp)) deallocate(V_interp)
if (allocated(c_eim_fast)) deallocate(c_eim_fast)
if (allocated(V_interp_sel)) deallocate(V_interp_sel)
call emulator_free(emu)

! Free physics arrays
if (allocated(k_ch)) deallocate(k_ch)
if (allocated(eta_val)) deallocate(eta_val)
if (allocated(E_thresh)) deallocate(E_thresh)
if (allocated(l_ch_arr)) deallocate(l_ch_arr)
if (allocated(ch_open_arr)) deallocate(ch_open_arr)
if (allocated(u_couple_arr)) deallocate(u_couple_arr)
if (allocated(b_full)) deallocate(b_full)
if (allocated(c_exact)) deallocate(c_exact)
if (allocated(S_column)) deallocate(S_column)
if (allocated(x_mesh)) deallocate(x_mesh)
if (allocated(r_mesh_arr)) deallocate(r_mesh_arr)
if (allocated(weights)) deallocate(weights)
if (allocated(lambda_arr)) deallocate(lambda_arr)
if (allocated(T_matrix)) deallocate(T_matrix)
if (allocated(atilde_J)) deallocate(atilde_J)
if (allocated(k_ch_all)) deallocate(k_ch_all)
if (allocated(eta_all)) deallocate(eta_all)
if (allocated(l_ch_all)) deallocate(l_ch_all)

! Free output file arrays
if (allocated(psi_exact_r)) deallocate(psi_exact_r)
if (allocated(psi_matrix_r)) deallocate(psi_matrix_r)
if (allocated(sigma_J_exact)) deallocate(sigma_J_exact)
if (allocated(sigma_J_matrix)) deallocate(sigma_J_matrix)
if (allocated(S11_J_exact)) deallocate(S11_J_exact)
if (allocated(S11_J_matrix)) deallocate(S11_J_matrix)
if (allocated(x_mesh_stored)) deallocate(x_mesh_stored)
if (allocated(S_matrix)) deallocate(S_matrix)

end program emulator_predict
