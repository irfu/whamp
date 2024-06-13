module comcout
   implicit none
!integer, parameter,private :: d2p=kind(1.0d0)
   integer, parameter, private :: d2p = 8

   complex(kind=d2p) :: X       ! frequency
   complex(kind=d2p) :: EFL(3)  ! E field
   complex(kind=d2p) :: BFL(3)  ! B field
   complex(kind=d2p) :: D       ! dispersion function
   complex(kind=d2p) :: DX      ! D derivative wrt frequency
   complex(kind=d2p) :: DZ      ! D derivative wrt Z
   complex(kind=d2p) :: DP      ! D derivative wrt P
   complex(kind=d2p) :: E(6, 4) ! epsilon and its derivatives x eps_x, z eps_z, p eps_p
   complex(kind=d2p) :: RI      ! complex refractive index

   real(kind=d2p) :: P
   real(kind=d2p) :: Z
   real(kind=d2p) :: VG(2)
   real(kind=d2p) :: SG(2)
   real(kind=d2p) :: ENE

end module comcout
