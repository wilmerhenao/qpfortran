!  This is a program that finds the solution to the QP problem
!  []
! [x,d,q,info] = qpspe!ial(G,varargin)
!
! Solves the QP
!
!    min     q(x)  = || G*x ||_2^2 = x'*(G'*G)*x
!    s.t.  sum(x)  = 1
!              x  >= 0
!
! The problem !orresponds to finding the smallest ve!tor
! (2-norm) in the !onvex hull of the !olumns of G
!
! Inputs:
!     G            -- (M x n double) matrix G, see problem above
!     varargin{1}  -- (int) maximum number of iterates allowed
!                     If not present, maxit = 100 is used
!     varargin{2}  -- (n x 1 double) ve!tor x0 with initial (FEASIBLE) iterate.
!                     If not present, (or requirements on x0 not met) a
!                     useable default x0 will be used
!
! Outputs:
!     x       -- Optimal point attaining optimal value
!     d = G*x -- Smallest ve!tor in the !onvex hull
!     q       -- Optimal value found = d'*d
!     info    -- Run data:
!                info(1) =
!                   0 = everything went well, q is optimal
!                   1 = maxit rea!hed and final x is feasible. so q
!                       might not be optimal, but it is better than q(x0)
!                   2 = something went wrong
!                info(2) = #iterations used

subroutine qpspecial(m, n, G, maxit)
implicit none
integer          m, n, maxit
double precision, dimension(m, n) :: G
double precision, dimension(n, n) :: invTrC

!  This is a program that finds the solution to the QP problem
!  []
! [x,d,q,info] = qpspecial(G,varargin)
!
! Solves the QP
!
!    min     q(x)  = || G*x ||_2^2 = x'*(G'*G)*x
!    s.t.  sum(x)  = 1
!              x  >= 0
!
! The problem !orresponds to finding the smallest vector
! (2-norm) in the convex hull of the columns of G
!
! Inputs:
!     G            -- (M x n double) matrix G, see problem above
!     varargin{1}  -- (int) maximum number of iterates allowed
!                     If not present, maxit = 100 is used
!     varargin{2}  -- (n x 1 double) vector x0 with initial (FEASIBLE) iterate.
!                     If not present, (or requirements on x0 not met) a
!                     useable default x0 will be used
!
! Outputs:
!     x       -- Optimal point attaining optimal value
!     d = G*x -- Smallest vector in the convex hull
!     q       -- Optimal value found = d'*d
!     info    -- Run data:
!                info(1) =
!                   0 = everything went well, q is optimal
!                   1 = maxit reached and final x is feasible. so q
!                       might not be optimal, but it is better than q(x0)
!                   2 = something went wrong
!                info(2) = #iterations used

integer          :: echo, info, k
double precision :: ptemp, eta, delta, mu0, tolmu, tolrs, kmu, nQ, krs, ap, ad, M, r2, rs, mu
double precision :: r5, r6, dy, ptemp, muaff
double precision, dimension(n)    :: x, y, z, zdx, KT, r1, r3, r4, r7, e, work, dx, dz, p
double precision, dimension(n, n) :: Q, QD, C, invTrC, invC
integer, dimension(n) :: ipiv, idx

! External procedures defined in lapack
external DGETRF
external DGETRI
echo = 0
! Check the dimensions

if (m*n .le. 0) then
   task = 'STOP:  Error in the dimensions of G'
   write(*,*) 'qpspecial is empty'
endif

! Create a vector of ones and use it as a starting point

do 1000 i = 1, n
   e(i) = 1d0
1000 continue
   
x = e
   
do 2100 i = 1, n
   idx(i) = (i - 1) * (n + 1) + 1
2100 continue
      
Q = matmul(transpose(G), G)

z = x
y = 0d0
eta = 0.9995d0
delta = 3d0
mu0 = dot_product(x, z) / n
tolmu = 1d-5
tolrs = 1d-5
kmu = tolmu * mu0
nQ = norminf(Q) + 2
krs = tolrs * nQ
ap = 0d0
ad = 0d0

do 2122 k = 1, maxit
   r1 = -matmul(Q,x) + e*y + z
   r2 = -1d0 + SUM(x)
   r3 = -x*z   ! double check this part
   rs = MAX(sum(abs(r1)), sum(abs(r2)))
   mu = -sum(r3)/n
   if(mu .lt. kmu) then
      if(rs .lt. krs) then
         goto 999
      end if
   end if
   zdx = z / x
   QD = Q
   QD(idx) = QD(idx) + zdx
   cholesky_sub(QD)
   C = QD
   invTrC = transpose(C)
   call DGETRF(n, n, invTrC, n, ipiv, info)
   if (0 /= info) then
      stop 'Matrix is numerically singular!'
   endif
   
   call DGETRF(n, n, invC, n, ipiv, info)
   if (0 /= info) then
      stop 'Matrix is numerically singular!'
   endif
   
   call DGETRI(n, invTrC, n, ipiv, work, n, info)
   if(0 /= info) then
      stop 'Matrix inversion failed!'
   endif

   call DGETRI(n, invC, n, ipiv, work, n, info)
   if(0 /= info) then
      stop 'Matrix inversion failed!'
   endif
   
   KT = matmul(invTrC, e)
   M = dot_product(KT, KT)
   
   r4 = r1 + r3 / x
   r5 = matmul(transpose(KT), matmul(invTrC, r4))
   r6 = r2 + r5
   dy = -r6 / M
   r7 = r4 + e * dy
   dx = matmul(invC, matmul(invTrC, r7))
   dz = (r3 - z * dx) / x
   p = -x / dx
   ptemp = 0d0
   do 2142 i = 1, n
      if (p(i) .gt. 0d0) then 
         min(p(i), ptemp)
      endif
2142  continue
   
   ap = min(ptemp, 1d0)
   p = -z / dz
   
do 2143 i = 1, n
      if (p(i) .gt. 0d0) then 
         min(p(i), ptemp)
      endif
2143  continue
   ad = min(ptemp, 1d0)
   
   muaff = matmul(transpose(x + ap * dx), (z + ad * dz)) / n
   sig = (muaff/mu)**delta
   r3 = r3 + sig * mu
   r3 = r3 - dx * dz
   r4 = r1 + r3 / x
   r5 = matmul(transpose(KT), matmul(invTrC, r4))
   r6 = r2 + r5
   dy = -r6/M
   r7 = r4 + matmul(e, dy)
   dx = matmul(invC, matmul(invTrC, r7))
   dz = (r3-z*dx)/x

   p = -x/dx
   do 2144 i = 1, n
      if (p(i) .gt. 0d0) then 
         min(p(i), ptemp)
      endif
2144  continue
   ap = min(ptemp, 1d0)
   
   p = -z / dz
   do 2145 i = 1, n
      if (p(i) .gt. 0d0) then 
         min(p(i), ptemp)
      endif
2145  continue
   ad = min(ptemp, 1d0)

   x=x+eta*ap*dx
   y=y+eta*ad*dy
   z=z+eta*ad*dz
   
2122 continue
x = max(x, 0d0)
x = x/sum(x)
d = matmul(G, x)
q = matmul(transpose(d),d)

999
end

! ----------------------- end of qpspecial ---------------------------------

function norminf(m, n, G)
integer m, n
double precision G(m, n), normvalue, normvaluetemp

normvalue = 0

do 1100 i = 1, m
   normvaluetemp = 0d0
   do 1200 j = 1, n
      normvaluetemp = normvaluetemp + abs(G(i, j))
   1200  continue
      normvalue = max(normvalue, normvaluetemp)
1100  continue
      normvalue
end function norminf
! ----------------------- end of norminf ------------------------------------

subroutine cholesky_sub(A,n)

implicit none

! formal vars
integer :: n      ! number of rows/cols in matrix
real    :: A(n,n) ! matrix to be decomposed

! local vars
integer :: j      ! iteration counter

! begin loop
chol: do j = 1,n

! perform diagonal component
A(j,j) = sqrt(A(j,j) - dot_product(A(j,1:j-1),A(j,1:j-1)))

! perform off-diagonal component
if (j < n) A(j+1:n,j) = (A(j+1:n,j) - matmul(A(j+1:n,1:j-1),A(j,1:j-1))) / &
&           A(j,j)

end do chol

end subroutine cholesky_sub

! ------------------------ end of cholesky_sub -----------------------------
