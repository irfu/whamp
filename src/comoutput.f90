module comoutput
   implicit none
!integer, parameter,private :: d2p=kind(1.0d0)
   integer, parameter, private :: d2p = 8
   integer  :: kperpSize, kparSize

   real(kind=d2p), allocatable :: kperpOUT(:) ! ???? should those be matrices or vectors?
   real(kind=d2p), allocatable :: kparOUT(:)

   complex(kind=d2p), allocatable :: fOUT(:, :)
   complex(kind=d2p), allocatable :: ExOUT(:, :)
   complex(kind=d2p), allocatable :: EyOUT(:, :)
   complex(kind=d2p), allocatable :: EzOUT(:, :)
   complex(kind=d2p), allocatable :: BxOUT(:, :)
   complex(kind=d2p), allocatable :: ByOUT(:, :)
   complex(kind=d2p), allocatable :: BzOUT(:, :)
   complex(kind=d2p), allocatable :: SxOUT(:, :)
   complex(kind=d2p), allocatable :: SyOUT(:, :)
   complex(kind=d2p), allocatable :: SzOUT(:, :)

   real(kind=d2p), allocatable :: EBOUT(:, :)
   real(kind=d2p), allocatable :: VGPOUT(:, :)
   real(kind=d2p), allocatable :: VGZOUT(:, :)
   real(kind=d2p), allocatable :: SGPOUT(:, :)
   real(kind=d2p), allocatable :: SGZOUT(:, :)
   real(kind=d2p), allocatable :: uOUT(:, :)

   integer, allocatable :: flagSolutionFoundOUT(:, :)
   integer, allocatable :: flagTooHeavilyDampedOUT(:, :)
   integer, allocatable :: flagNoConvergenceOUT(:, :)

end module comoutput
