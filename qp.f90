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
double precision G(m, n)


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

integer echo, idx(n)
double precision e(n), x(n)

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
tolmu = 1e-5
tolrs = 1e-5
kmu = tolmu * mu0
nQ = norm(Q, inf) + 2
