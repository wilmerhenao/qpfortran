program driver
  implicit none
  double precision, dimension(4, 3) :: G
  double precision, dimension(3) :: myx
  integer function qpspecial(m, n, G, maxit)


G(1, 1) = 1d0 
G(1, 2) = 2d0
G(1, 3) = 3d0
G(2, 1) = 4d0
G(2, 2) = 3.2d0
G(2, 3) = 21d0
G(3, 1) = 2.21d0
G(3, 2) = 0d0
G(3, 3) = 1.122321d0;
G(4, 1) = 2.2d0
G(4, 2) = -3.1d0
G(4, 3) = 0.5d0

myx = qpspecial(4, 3, G, 100)

end program driver
