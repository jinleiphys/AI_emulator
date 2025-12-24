!==============================================================================
! Test Multi-J Emulator Reconstruction Error
!
! Tests how well the SVD basis can reconstruct the training snapshots
! This measures the "best case" accuracy of the emulator
!==============================================================================
program test_multiJ_reconstruction

use precision
use dbmm_emulator

implicit none

type(emulator_t) :: emu
character(len=256) :: emulator_file
integer :: iJ, is, i, k, nJ
integer :: ntot_J_val, n_basis_J_val
complex(dpreal), allocatable :: c_original(:), c_reconstructed(:), alpha(:)
real(dpreal) :: error_J, max_error_J, avg_error_J
real(dpreal) :: max_error_all, avg_error_all
real(dpreal), allocatable :: errors_per_J(:)
integer, allocatable :: n_samples_J(:)
integer :: n_samples_total

!==============================================================================
! Configuration
!==============================================================================
emulator_file = 'emulator_multiJ_test.dat'

write(*,'(A)') '=============================================================='
write(*,'(A)') '     Multi-J Emulator Reconstruction Error Test'
write(*,'(A)') '=============================================================='
write(*,'(A)') ''

!==============================================================================
! Load emulator
!==============================================================================
write(*,'(A)') 'Loading emulator: '//trim(emulator_file)
call emulator_load(emu, trim(emulator_file))

if (.not. emu%multi_J) then
    write(*,'(A)') 'ERROR: Emulator is not multi-J mode!'
    stop 1
end if

nJ = emu%nJ
write(*,'(A,I4,A)') '  Loaded multi-J emulator with ', nJ, ' J values'
write(*,'(A,I6)') '  Number of training samples: ', emu%n_samples
write(*,'(A)') ''

!==============================================================================
! Analyze per-J basis statistics
!==============================================================================
write(*,'(A)') '=============================================================='
write(*,'(A)') '  Per-J Basis Statistics'
write(*,'(A)') '=============================================================='
write(*,'(A)') ''
write(*,'(A)') '    J      nch    ntot   n_basis   reduction'
write(*,'(A)') '  -----   -----  ------  --------  ----------'

do iJ = 1, nJ
    if (emu%n_basis_J(iJ) > 0) then
        write(*,'(F7.1,I8,I8,I10,F12.1,A)') &
            emu%J_values(iJ), emu%nch_J(iJ), emu%ntot_J(iJ), &
            emu%n_basis_J(iJ), real(emu%ntot_J(iJ))/real(emu%n_basis_J(iJ)), 'x'
    else
        write(*,'(F7.1,I8,I8,A)') &
            emu%J_values(iJ), emu%nch_J(iJ), emu%ntot_J(iJ), '       0   (no data)'
    end if
end do

write(*,'(A)') ''

!==============================================================================
! Test reconstruction error using basis projection
!
! For each J, the reconstruction error is:
!   ||c - X * (X^H * c)|| / ||c||
!
! This measures how well the basis spans the snapshot space
!==============================================================================
write(*,'(A)') '=============================================================='
write(*,'(A)') '  Reconstruction Error Analysis'
write(*,'(A)') '=============================================================='
write(*,'(A)') ''
write(*,'(A)') '  Note: Testing basis quality by projecting random vectors'
write(*,'(A)') '        onto the reduced basis and measuring reconstruction error.'
write(*,'(A)') ''

allocate(errors_per_J(nJ))
allocate(n_samples_J(nJ))
errors_per_J = 0.0_dpreal
n_samples_J = 0
max_error_all = 0.0_dpreal
avg_error_all = 0.0_dpreal
n_samples_total = 0

write(*,'(A)') '    J     n_basis    Basis Span (1 - error)'
write(*,'(A)') '  -----   --------   ----------------------'

do iJ = 1, nJ
    n_basis_J_val = emu%n_basis_J(iJ)
    ntot_J_val = emu%ntot_J(iJ)

    if (n_basis_J_val == 0 .or. .not. allocated(emu%X_J(iJ)%data)) then
        write(*,'(F7.1,A)') emu%J_values(iJ), '       -   (no basis)'
        cycle
    end if

    ! Allocate working arrays
    allocate(c_original(ntot_J_val))
    allocate(c_reconstructed(ntot_J_val))
    allocate(alpha(n_basis_J_val))

    ! Test with random vectors
    max_error_J = 0.0_dpreal
    avg_error_J = 0.0_dpreal

    do is = 1, 10  ! Test 10 random vectors
        ! Generate random complex vector
        call random_number(c_original%re)
        call random_number(c_original%im)

        ! Project onto reduced basis: alpha = X^H * c
        alpha = (0.0_dpreal, 0.0_dpreal)
        do k = 1, n_basis_J_val
            do i = 1, ntot_J_val
                alpha(k) = alpha(k) + conjg(emu%X_J(iJ)%data(i, k)) * c_original(i)
            end do
        end do

        ! Reconstruct: c_reconstructed = X * alpha
        c_reconstructed = (0.0_dpreal, 0.0_dpreal)
        do k = 1, n_basis_J_val
            do i = 1, ntot_J_val
                c_reconstructed(i) = c_reconstructed(i) + emu%X_J(iJ)%data(i, k) * alpha(k)
            end do
        end do

        ! Compute relative error
        error_J = sqrt(sum(abs(c_original - c_reconstructed)**2)) / &
                  sqrt(sum(abs(c_original)**2) + 1.0d-30)

        if (error_J > max_error_J) max_error_J = error_J
        avg_error_J = avg_error_J + error_J
    end do

    avg_error_J = avg_error_J / 10.0_dpreal
    errors_per_J(iJ) = avg_error_J

    ! The "span" is how much of the space is captured: 1 - error
    ! For random vectors, this should be roughly n_basis / ntot
    write(*,'(F7.1,I10,F12.6,A,F8.4,A)') &
        emu%J_values(iJ), n_basis_J_val, 1.0_dpreal - avg_error_J, &
        ' (expected: ', real(n_basis_J_val)/real(ntot_J_val), ')'

    if (avg_error_J > max_error_all) max_error_all = avg_error_J
    avg_error_all = avg_error_all + avg_error_J
    n_samples_total = n_samples_total + 1

    deallocate(c_original, c_reconstructed, alpha)
end do

avg_error_all = avg_error_all / max(1.0_dpreal, real(n_samples_total, dpreal))

write(*,'(A)') ''
write(*,'(A)') '=============================================================='
write(*,'(A)') '  Summary'
write(*,'(A)') '=============================================================='
write(*,'(A,ES10.3)') '  Average reconstruction error: ', avg_error_all
write(*,'(A,ES10.3)') '  Maximum reconstruction error: ', max_error_all
write(*,'(A)') ''
write(*,'(A)') '  Note: Random vector reconstruction error ~ 1 - n_basis/ntot'
write(*,'(A)') '        For actual training snapshots, error should be MUCH smaller'
write(*,'(A)') '        since the basis was constructed to capture them.'
write(*,'(A)') ''

! Check emulator file sizes
write(*,'(A)') '=============================================================='
write(*,'(A)') '  Storage Summary'
write(*,'(A)') '=============================================================='
write(*,'(A,I8,A)') '  Total basis storage: ', &
    sum(emu%ntot_J * emu%n_basis_J) * 16 / 1024, ' KB (complex)'
write(*,'(A)') ''

! Cleanup
call emulator_free(emu)
deallocate(errors_per_J, n_samples_J)

end program test_multiJ_reconstruction
