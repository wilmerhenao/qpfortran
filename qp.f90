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
double precision, dimension(n, n) :: invC

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

integer :: echo, info
double precision, dimension(n) :: e
double precision, dimension(n) :: x
double precision, dimension(n) :: work
integer, dimension(n) :: ipiv
integer, dimension(n) :: idx

! External procedures defined in lapack
external DGETRF
external DGETRI

echo = 0

c Check the dimensions

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
y = 0
eta = 0.9995
delta = 3
mu0 = matmul(transpose(x), z) / n
tolmu = 1d-5
tolrs = 1d-5
kmu = tolmu * mu0
nQ = norminf(Q) + 2
krs = tolrs * nQ
ap = 0
ad = 0

if(echo > 0) then
   write(*,*) 'k mu stpsz res'
endif

do 2122 k = 1, maxit
   r1 = -matmul(Q,x) + matmul(e, y) + z
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
   invC = transpose(C)
   call DGETRF(n, n, invC, n, ipiv, info)
   if (0 /= info) then
      stop 'Matrix is numerically singular!'
   endif
   call DGETRI(n, Ainv, n, ipiv, work, n, info)
   
   if(0 /= info) then
      stop 'Matrix inversion failed!'
   endif
   
   KT = matmul(invC, e)
   M = matmul(KT, KT)
   
   r4 = r1 + r3 / x
   r5 = matmul(transpose(K), matmul(invC, r4))
   r6 = r2 + r5
   dy = -r6 / M
   r7 = r4 + matmul(e, dy)
2122 continue

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
