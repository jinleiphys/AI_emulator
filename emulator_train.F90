!==============================================================================
! Emulator Training Program for CDCC Calculations
!
! Usage: ./emulator_train <config_file>
! Example: ./emulator_train emulator_train_dNi.in
!
! This program:
! 1. Reads training configuration (parameter ranges, n_samples, etc.)
! 2. Initializes the physics system from base_input
! 3. For each LHS sample: modifies potentials, runs CDCC, collects coefficients
! 4. Trains the emulator (SVD or Gram-Schmidt)
! 5. Saves the trained emulator
!==============================================================================
program emulator_train

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

type(train_config_t) :: config
type(emulator_t) :: emu
integer :: ierr, i, j, is, njtot, parity, ich

! Training data
real(dpreal), allocatable :: param_samples(:,:)
real(dpreal), allocatable :: param_min(:), param_max(:)
real(dpreal), allocatable :: base_params(:)
character(len=32), allocatable :: param_names(:)

! CDCC variables
type(nch_tilde) :: atilde
type(nch_tilde), allocatable :: atilde_J(:)  ! Multi-J: atilde for each J
real(dpreal) :: mu_val, ecm_val, jtot
integer :: ngs, ntot_max
real(dpreal), allocatable :: k_ch(:), eta_val(:), E_thresh(:)
integer, allocatable :: l_ch_arr(:), ch_open_arr(:)
complex(dpreal), allocatable :: u_couple_arr(:,:,:)
complex(dpreal), allocatable :: S_column(:), c_vec(:)

! Multi-J variables
integer :: nJ, iJ, nch_max_all
real(dpreal), allocatable :: J_values(:)
integer, allocatable :: nch_J(:)

! Random J sampling variables
integer :: n_J_per_sample, iJ_sample
integer, allocatable :: J_sample_indices(:)

! Timing
real(dpreal) :: t_start, t_end, t_sample_start, t_sample_end
real(dpreal) :: time_total, time_per_sample
integer :: n_collected

!==============================================================================
! Step 0: Read training configuration from command line
!==============================================================================
write(*,'(A)') '=============================================================='
write(*,'(A)') '        DBMM Emulator Training Program'
write(*,'(A)') '=============================================================='
write(*,'(A)') ''

! Read config directly from stdin
call read_train_config_stdin(config, ierr)
if (ierr /= 0) then
    write(*,'(A)') 'ERROR: Failed to read training configuration'
    stop 1
end if

!==============================================================================
! Step 1: Initialize physics system from base_input
!==============================================================================
write(*,'(A)') ''
write(*,'(A)') 'Step 1: Initializing physics system...'
write(*,'(A)') '  Base input: '//trim(config%base_input)

! Redirect stdin to read from base_input
open(unit=5, file=trim(config%base_input), status='old', iostat=ierr)
if (ierr /= 0) then
    write(*,'(A)') 'ERROR: Cannot open base input: '//trim(config%base_input)
    stop 1
end if

call initialize()

! Setup bins and wavefunctions (need to continue reading from unit 5)
if (method == 1) then
    call setup_bins()
else
    call tho_wavefunction()
end if

! Close input file after all reading is done
close(5)

call alpha_3b_new()
call phibx_bound()
call phibx_bin()

write(*,'(A)') '  Physics system initialized'
write(*,'(A,I5)') '  Number of 3-body channels: ', a3b%nchmax
write(*,'(A,I5)') '  nlag: ', nlag

!==============================================================================
! Step 2: Setup parameter sampling
!==============================================================================
write(*,'(A)') ''
write(*,'(A)') 'Step 2: Setting up parameter sampling...'

allocate(param_names(config%n_params))
allocate(param_min(config%n_params))
allocate(param_max(config%n_params))
allocate(base_params(config%n_params))
allocate(param_samples(config%n_params, config%n_samples))

do i = 1, config%n_params
    param_names(i) = config%params(i)%name
    param_min(i) = config%params(i)%min_val
    param_max(i) = config%params(i)%max_val
    base_params(i) = get_potential_param(param_names(i), ierr)
    if (ierr /= 0) then
        write(*,'(A,A)') '  WARNING: Could not get base value for: ', trim(param_names(i))
        base_params(i) = (param_min(i) + param_max(i)) / 2.0_dpreal
    end if
    write(*,'(A,A,A,F10.3,A,F10.3,A,F10.3,A)') '  ', trim(param_names(i)), &
        ': base=', base_params(i), ', range=[', param_min(i), ', ', param_max(i), ']'
end do

! Generate LHS samples
call lhs_sample(config%n_samples, config%n_params, param_min, param_max, param_samples)
write(*,'(A,I5,A)') '  Generated ', config%n_samples, ' LHS samples'

!==============================================================================
! Step 3: Determine maximum system size and initialize emulator
!==============================================================================
write(*,'(A)') ''
write(*,'(A)') 'Step 3: Determining system dimensions...'

! Calculate ecm and mu
ecm_val = elab * masst / (massp + masst)
mu_val = amu * (masst * massp) / (massp + masst)

ecm_val = ecm_val + epsilon_bin(1,1)

! Initialize emulator
if (config%n_basis == 0) then
    ! Auto-determine basis size (will be set by GS method)
    config%n_basis = min(config%n_samples, 50)
end if

if (config%multi_J) then
    !-----------------------------------------------------------------------
    ! Multi-J mode: scan all J values from 0 to jmax
    !-----------------------------------------------------------------------
    write(*,'(A)') '  Multi-J mode enabled'

    ! Determine number of J values (J = 0, 1, 2, ..., jmax)
    nJ = jmax + 1
    allocate(J_values(nJ))
    allocate(nch_J(nJ))
    allocate(atilde_J(nJ))

    nch_max_all = 0

    write(*,'(A,I4,A,I4,A)') '  Scanning J = 0 to ', int(jmax), ', total ', nJ, ' values'

    do iJ = 1, nJ
        jtot = real(iJ - 1, dpreal)  ! J = 0, 1, 2, ..., jmax
        J_values(iJ) = jtot
        ! For zero spin: parity = (-1)^J (even J -> +1, odd J -> -1)
        parity = 1 - 2 * mod(iJ - 1, 2)

        call alphatildeindex(jtot, parity, atilde_J(iJ))
        nch_J(iJ) = atilde_J(iJ)%nchmax

        ! Compute nch_open for each J
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

        if (nch_J(iJ) > nch_max_all) nch_max_all = nch_J(iJ)

        if (mod(iJ, 5) == 1 .or. iJ == nJ) then
            write(*,'(A,F5.1,A,I4,A,I4)') '    J=', jtot, ': nch=', nch_J(iJ), &
                ', nch_open=', atilde_J(iJ)%nch_open
        end if
    end do

    write(*,'(A,I5)') '  Max channels across all J: ', nch_max_all
    write(*,'(A,I5)') '  Max ntot: ', nlag * nch_max_all

    ! Initialize emulator (per-J basis)
    write(*,'(A)') '  Using PER-J separate bases'
    call emulator_init_multiJ(emu, nlag, nJ, J_values, nch_J, config%n_basis, config%n_params)
    call emulator_allocate_training_multiJ(emu, config%n_samples)

    ! Use first J's atilde for compatibility with later code
    atilde = atilde_J(1)

else
    !-----------------------------------------------------------------------
    ! Single-J mode (backward compatible)
    !-----------------------------------------------------------------------
    write(*,'(A)') '  Single-J mode (first J only)'

    ! Find maximum ntot by checking first J
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
        j = 0
        do ich = 1, atilde%nchmax
            if ((ecm_val - epsilon_bin(atilde%n(ich), a3b%alpha2b(atilde%alpha3b(ich)))) > -0.00001_dpreal) then
                j = j + 1
                atilde%ch_open(j) = ich
            end if
        end do
    end if

    ntot_max = nlag * atilde%nchmax

    write(*,'(A,I5)') '  First J nchmax: ', atilde%nchmax
    write(*,'(A,I5)') '  ntot for first J: ', ntot_max

    call emulator_init(emu, nlag, atilde%nchmax, config%n_basis, config%n_params)
    call emulator_allocate_training(emu, config%n_samples)

    ! Also allocate EIM training storage for potential snapshots
    call eim_allocate_training(emu, config%n_samples)
end if

write(*,'(A,I5)') '  Emulator initialized with n_basis: ', config%n_basis

!==============================================================================
! Step 4: Training loop - collect snapshots
!==============================================================================
write(*,'(A)') ''
write(*,'(A)') 'Step 4: Collecting training snapshots...'
write(*,'(A)') ''

call cpu_time(t_start)
n_collected = 0

if (config%multi_J) then
    !-----------------------------------------------------------------------
    ! Multi-J training loop
    !-----------------------------------------------------------------------

    ! Determine J sampling strategy
    if (config%random_J) then
        ! Random J sampling: select n_J_samples J values per parameter set
        if (config%n_J_samples <= 0 .or. config%n_J_samples > nJ) then
            n_J_per_sample = nJ  ! Default to all J values
        else
            n_J_per_sample = config%n_J_samples
        end if
        write(*,'(A,I4,A,I4,A)') '  Random J sampling: ', n_J_per_sample, &
            ' J values per sample x ', config%n_samples, ' samples'
        allocate(J_sample_indices(n_J_per_sample))
    else
        ! Full J scan: compute all J values for each parameter set
        n_J_per_sample = nJ
        write(*,'(A,I4,A,I4,A)') '  Full J scan: ', nJ, ' J values x ', config%n_samples, ' samples'
    end if

    ! Allocate working arrays using max dimensions
    allocate(k_ch(nch_max_all))
    allocate(eta_val(nch_max_all))
    allocate(E_thresh(nch_max_all))
    allocate(l_ch_arr(nch_max_all))
    allocate(ch_open_arr(nch_max_all))
    allocate(u_couple_arr(nlag, nch_max_all, nch_max_all))
    allocate(S_column(nch_max_all))
    allocate(c_vec(nlag * nch_max_all))

    do is = 1, config%n_samples
        call cpu_time(t_sample_start)

        ! Set potential parameters for this sample
        do i = 1, config%n_params
            call set_potential_param(param_names(i), param_samples(i, is), ierr)
        end do

        ! Recompute potentials with new parameters
        call potr('t', 1, zb*zt, 0.0_dpreal)
        UbA = V
        call potr('x', 1, zx*zt, 0.0_dpreal)
        UxA = V

        ! Select J values for this sample
        if (config%random_J) then
            ! Randomly select n_J_per_sample distinct J indices
            call random_select_J(nJ, n_J_per_sample, J_sample_indices)
        end if

        ! Loop over selected J values
        do iJ_sample = 1, n_J_per_sample
            if (config%random_J) then
                iJ = J_sample_indices(iJ_sample)
            else
                iJ = iJ_sample
            end if
            ! Clear potential cache for fresh computation
            call cleanup_potential_cache()

            ! Setup channel information for this J
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

            ! Solve DBMM for incoming channel 1
            c_vec = (0.0_dpreal, 0.0_dpreal)
            call dbmm_coupled_channel_cc(atilde_J(iJ)%nchmax, atilde_J(iJ)%nch_open, &
                                          ch_open_arr(1:atilde_J(iJ)%nch_open), &
                                          nlag, rmax, l_ch_arr(1:atilde_J(iJ)%nchmax), &
                                          k_ch(1:atilde_J(iJ)%nchmax), &
                                          eta_val(1:atilde_J(iJ)%nchmax), &
                                          u_couple_arr(1:nlag, 1:atilde_J(iJ)%nchmax, 1:atilde_J(iJ)%nchmax), &
                                          mu_val, ecm_val, E_thresh(1:atilde_J(iJ)%nchmax), &
                                          1, S_column(1:atilde_J(iJ)%nch_open), &
                                          c_out=c_vec(1:nlag*atilde_J(iJ)%nchmax), &
                                          solver_type=solver_type)

            ! Store snapshot for this J (per-J basis)
            call emulator_store_snapshot_J(emu, iJ, is, c_vec(1:nlag*atilde_J(iJ)%nchmax), &
                                           param_samples(:, is))
            n_collected = n_collected + 1
        end do

        call cpu_time(t_sample_end)

        ! Progress output
        if (mod(is, 5) == 0 .or. is == 1 .or. is == config%n_samples) then
            write(*,'(A,I4,A,I4,A,I5,A,F8.1,A)') '  Sample ', is, '/', config%n_samples, &
                ' (', n_J_per_sample, ' J values, ', (t_sample_end - t_sample_start)*1000.0_dpreal, ' ms)'
        end if
    end do

    ! Clean up random J sampling arrays
    if (config%random_J .and. allocated(J_sample_indices)) then
        deallocate(J_sample_indices)
    end if

else
    !-----------------------------------------------------------------------
    ! Single-J training loop (backward compatible)
    !-----------------------------------------------------------------------
    ! Allocate working arrays
    allocate(k_ch(atilde%nchmax))
    allocate(eta_val(atilde%nchmax))
    allocate(E_thresh(atilde%nchmax))
    allocate(l_ch_arr(atilde%nchmax))
    allocate(ch_open_arr(atilde%nch_open))
    allocate(u_couple_arr(nlag, atilde%nchmax, atilde%nchmax))
    allocate(S_column(atilde%nch_open))
    allocate(c_vec(nlag * atilde%nchmax))

    do is = 1, config%n_samples
        call cpu_time(t_sample_start)

        ! Set potential parameters for this sample
        do i = 1, config%n_params
            call set_potential_param(param_names(i), param_samples(i, is), ierr)
        end do

        ! Recompute potentials with new parameters
        call potr('t', 1, zb*zt, 0.0_dpreal)
        UbA = V
        call potr('x', 1, zx*zt, 0.0_dpreal)
        UxA = V

        ! CRITICAL: Clear potential cache so u_alpha_tilde uses new potentials
        call cleanup_potential_cache()

        ! Setup channel information
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

        do ich = 1, atilde%nch_open
            ch_open_arr(ich) = ich
        end do

        ! Compute coupling potential (DBMM handles kinetic/boundary internally)
        call u_alpha_tilde(u_couple_arr, atilde)

        ! Store potential snapshot for EIM training
        call eim_store_potential(emu, is, u_couple_arr)

        ! Solve DBMM for incoming channel 1
        call dbmm_coupled_channel_cc(atilde%nchmax, atilde%nch_open, ch_open_arr, &
                                      nlag, rmax, l_ch_arr, k_ch, eta_val, &
                                      u_couple_arr, mu_val, ecm_val, E_thresh, &
                                      1, S_column, c_out=c_vec, solver_type=solver_type)

        ! Store snapshot
        call emulator_store_snapshot(emu, is, c_vec, param_samples(:, is))
        n_collected = n_collected + 1

        call cpu_time(t_sample_end)

        ! Progress output
        if (mod(is, 10) == 0 .or. is == 1 .or. is == config%n_samples) then
            write(*,'(A,I4,A,I4,A,F8.3,A)') '  Sample ', is, '/', config%n_samples, &
                ' completed (', (t_sample_end - t_sample_start)*1000.0_dpreal, ' ms)'
        end if
    end do
end if

call cpu_time(t_end)
time_total = t_end - t_start
time_per_sample = time_total / real(config%n_samples, dpreal)

write(*,'(A)') ''
write(*,'(A,I5,A)') '  Collected ', n_collected, ' snapshots'
write(*,'(A,F10.3,A)') '  Total training data time: ', time_total, ' s'
write(*,'(A,F10.3,A)') '  Time per sample: ', time_per_sample * 1000.0_dpreal, ' ms'

! Restore original parameters
do i = 1, config%n_params
    call set_potential_param(param_names(i), base_params(i), ierr)
end do

!==============================================================================
! Step 5: Finalize training (SVD or Gram-Schmidt)
!==============================================================================
write(*,'(A)') ''
write(*,'(A)') 'Step 5: Training emulator...'

call cpu_time(t_start)

if (config%multi_J) then
    ! Multi-J finalization: SVD for each J independently
    write(*,'(A,ES10.2)') '  Method: Per-J SVD (tol=', config%tol, ')'
    call emulator_finalize_training_multiJ(emu, config%tol)
else
    if (trim(config%method) == 'gs') then
        write(*,'(A,ES10.2)') '  Method: Gram-Schmidt (auto basis, tol=', config%tol, ')'
        call emulator_finalize_training_gs(emu, config%tol)
    else
        write(*,'(A,I5)') '  Method: SVD (n_basis=', config%n_basis, ')'
        call emulator_finalize_training(emu)
    end if
end if

call cpu_time(t_end)

write(*,'(A,F10.3,A)') '  Training time: ', (t_end - t_start) * 1000.0_dpreal, ' ms'
if (config%multi_J) then
    write(*,'(A,I5,A)') '  Trained ', nJ, ' J-specific bases'
else
    write(*,'(A,I5)') '  Final basis size: ', emu%n_basis
end if

!==============================================================================
! Step 5b: Finalize EIM training for potential (single-J only for now)
!==============================================================================
if (.not. config%multi_J) then
    write(*,'(A)') ''
    write(*,'(A)') 'Step 5b: Training EIM for potential...'

    call cpu_time(t_start)

    ! Use same number of basis functions as wavefunction emulator, or auto-determine
    call eim_finalize_training(emu, emu%n_basis, config%tol)

    call cpu_time(t_end)
    write(*,'(A,F10.3,A)') '  EIM training time: ', (t_end - t_start) * 1000.0_dpreal, ' ms'
else
    write(*,'(A)') ''
    write(*,'(A)') 'Step 5b: EIM training skipped (multi-J mode)'
    ! TODO: Implement per-J EIM training
end if

!==============================================================================
! Step 6: Save emulator
!==============================================================================
write(*,'(A)') ''
write(*,'(A)') 'Step 6: Saving emulator...'
write(*,'(A)') '  Output file: '//trim(config%output_file)

call emulator_save(emu, trim(config%output_file))

write(*,'(A)') '  Emulator saved successfully'

!==============================================================================
! Summary
!==============================================================================
write(*,'(A)') ''
write(*,'(A)') '=============================================================='
write(*,'(A)') '                    TRAINING COMPLETE'
write(*,'(A)') '=============================================================='
write(*,'(A,I5)')     '  Training samples:     ', config%n_samples
write(*,'(A,I5)')     '  Parameters:           ', config%n_params
write(*,'(A,I5)')     '  Wavefunction basis:   ', emu%n_basis
write(*,'(A,I5)')     '  Original dimension:   ', emu%ntot
write(*,'(A,F10.1,A)')'  Dimension reduction:  ', real(emu%ntot)/real(emu%n_basis), 'x'
write(*,'(A,I5)')     '  EIM potential basis:  ', emu%n_pot_basis
write(*,'(A,I5)')     '  EIM interp points:    ', emu%n_pot_interp
write(*,'(A,A)')      '  Emulator file:        ', trim(config%output_file)
write(*,'(A)') '=============================================================='

! Cleanup
call emulator_free(emu)
deallocate(param_names, param_min, param_max, base_params, param_samples)
deallocate(k_ch, eta_val, E_thresh, l_ch_arr, ch_open_arr)
deallocate(u_couple_arr, S_column, c_vec)

contains

!==============================================================================
! Randomly select k distinct indices from 1 to n (Fisher-Yates shuffle variant)
!==============================================================================
subroutine random_select_J(n, k, indices)
    integer, intent(in) :: n, k
    integer, intent(out) :: indices(k)

    integer :: i, j, temp
    integer, allocatable :: pool(:)
    real :: r

    ! Initialize random seed (once)
    call random_seed()

    ! Create pool of all indices
    allocate(pool(n))
    do i = 1, n
        pool(i) = i
    end do

    ! Select k random indices using partial Fisher-Yates shuffle
    do i = 1, k
        ! Generate random index from i to n
        call random_number(r)
        j = i + int(r * real(n - i + 1))
        if (j > n) j = n

        ! Swap pool(i) and pool(j)
        temp = pool(i)
        pool(i) = pool(j)
        pool(j) = temp

        ! Store selected index
        indices(i) = pool(i)
    end do

    deallocate(pool)

end subroutine random_select_J

end program emulator_train
