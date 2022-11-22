program calculate_sph_nodes

use msis_init, only : msisinit

implicit none

integer :: nrec
integer :: num_sou
integer :: num_sta
integer :: num_doy
integer :: num_sec
! to deal with NRLMSIS2.0 input/output
integer :: mass
real(4) :: ap_msis(7), apd
real(4) :: d(9),t(2)
! to deal with HWM14
real(4) :: ap_hwm(2), w(2), aph
real(4), external :: pershift
! shared between NRMSLSIS2.0 and HWM14
integer :: iyd
real(4) :: sec, alt, glat, glon, stl, f107a, f107


integer, dimension(:), allocatable :: doy_vec, sec_vec

! dummy variables
integer :: idoy, isou, ista, ilin, num_lin, isec
integer :: OpenStatus

character(61) :: nodes_in, nodes_out
character(55) :: param_in
character(100) :: texty_text

!Initialize NRMLSIS2.0 model
call msisinit(parmpath='./',parmfile='msis20.parm')


write(*,*) "Reading ", "../input/secs.txt"
open (10, FILE="../input/secs.txt", STATUS="OLD", ACTION="read",        &
      POSITION="rewind", IOSTAT=OpenStatus)
if (OpenStatus > 0) STOP "*** Cannot open file ***"
num_sec = Count_Lines(10)
allocate(sec_vec(num_sec))
do ilin = 1, num_sec
    read(10, '(I3)') sec_vec(ilin)
enddo
close(10)

write(*,*) "Reading ", "../input/doys.txt"
open (10, FILE="../input/doys.txt", STATUS="OLD", ACTION="read",        &
      POSITION="rewind", IOSTAT=OpenStatus)
if (OpenStatus > 0) STOP "*** Cannot open file ***"
num_doy = Count_Lines(10)
allocate(doy_vec(num_doy))
do ilin = 1, num_doy
    read(10, '(I3)') doy_vec(ilin)
enddo
close(10)

write(*,*) "Reading ", "../input/sources.txt"
open (10, FILE="../input/sources.txt", STATUS="OLD", ACTION="read",        &
      POSITION="rewind", IOSTAT=OpenStatus)
if (OpenStatus > 0) STOP "*** Cannot open file ***"
num_sou = Count_Lines(10)
close(10)

write(*,*) "Reading ", "../input/stations.txt"
open (10, FILE="../input/stations.txt", STATUS="OLD", ACTION="read",        &
      POSITION="rewind", IOSTAT=OpenStatus)
if (OpenStatus > 0) STOP "*** Cannot open file ***"
num_sta = Count_Lines(10)
close(10)

! Input/output path+file name of each nodes file
100 format (A22, I5.5, A1, I3.3, A1, I5.5, A1, I4.4, A4)

do isec = 1, size(sec_vec)
    do idoy = 1, size(doy_vec)
        do isou = 1, num_sou
            do ista = 1, num_sta
                write(nodes_in, 100) "../output/nodes/nodes_", sec_vec(isec), "_", &
                                                            doy_vec(idoy), "_",  &
                                                            isou, "_",           &
                                                            ista, ".txt"
                write(*,*) "Reading ", trim(nodes_in)
                open (10, FILE=trim(nodes_in), STATUS="OLD", ACTION="read",        &
                    POSITION="rewind", IOSTAT=OpenStatus)
                if (OpenStatus > 0) STOP "*** Cannot open file ***"

                nrec = Count_Lines(10)

                write(nodes_out, 100) "../output/nodes/descr_", sec_vec(isec), "_", &
                                                                doy_vec(idoy),     &
                                                            "_", isou,          &
                                                            "_", ista, ".txt"

                ! Open file to save profiles with name created above
                open (20, FILE=TRIM(nodes_out), STATUS="REPLACE",              &
                    ACTION="write", POSITION="rewind", IOSTAT=OpenStatus)
                if (OpenStatus > 0) STOP "*** Cannot open file ***"
                write(*,*) "    Writing to ", trim(nodes_out)

                do ilin = 1, nrec
                    read(10,*) iyd,sec,alt,glat,glon,stl,f107a,f107,apd,aph
                    ap_msis(1) = apd ! daily ap
                    ap_hwm(2) = aph ! 3h ap indesx for 3 hrs before current time (fix this)
                    ! from checkwhm14.f90, not really used but maybe...
                    stl = pershift(sec + glon/15.0, (/0.0, 24.0/) )
                    call gtd8d(iyd,sec,alt,glat,glon,stl,f107a,f107,ap_msis,mass,d,t)
                    call hwm14(iyd,sec,alt,glat,glon,stl,f107a,f107,ap_hwm,w)
                    write(20,'(2i7,3f7.1,e13.4,f8.2,2f9.3)')  &
                    iyd,int(sec),alt,glat,glon,d(6),t(2),w(1),w(2)
                enddo
                close(20) ! nodes output file
                write(*,*) "    Done writing."
                close(10) ! nodes input file
                write(*,*) "Done reading."
            enddo
        enddo
    enddo
enddo


deallocate(doy_vec)

contains
    function Count_Lines(FileNum)
        integer, INTENT(IN) :: FileNum
        integer :: Count_Lines
        integer :: InputStatus2
        integer :: num_lin2
        num_lin2 = 0
        do
            read (FileNum, *, IOSTAT=InputStatus2)
            if (InputStatus2 < 0) EXIT ! End of file
            num_lin2 = num_lin2 + 1
        enddo
        Count_Lines = num_lin2
        rewind (FileNum)
    end function Count_Lines

end program calculate_sph_nodes

!******************************************************************************
!
!PERSHIFT
!JOHN EMMERT   9/12/03
!TRANSLATED TO FORTRAN-90 10/4/06. FORTRAN VERSION ONLY ALLOWS SCALAR INPUTS
!SHIFTS INPUT VALUES INTO A SPECIFIED PERIODIC INTERVAL
!
!CALLING SEQUENCE:   Result = PERSHIFT(x, range)
!
!ARGUMENTS
!      x:        The value to be shifted
!      perint:   2-element vector containing the start and end values
!                of the desired periodic interval.  The periodicity is
!                determined by the span of the range.
!
!ROUTINES USED THAT ARE NOT IN THE STANDARD FORTRAN-90 LIBRARY
!      None

function pershift(x, perint)

  real(4), parameter :: tol=1e-4
  real(4)            :: x, perint(0:1)
  real(4)            :: a, span, offset, offset1, pershift

  pershift = x
  a = perint(0)
  span = perint(1) - perint(0)
  if (span .ne. 0) then
    offset = x-a
    offset1 = mod(offset,span)
    if (abs(offset1) .lt. tol) offset1 = 0
  endif
  pershift = a + offset1
  if ((offset .lt. 0) .and. (offset1 .ne. 0)) pershift = pershift + span

  return

end function pershift
