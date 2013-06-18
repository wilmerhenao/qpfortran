c  This is a program that finds the solution to the QP problem
c  []
c [x,d,q,info] = qpspecial(G,varargin)
c
c Solves the QP
c
c    min     q(x)  = || G*x ||_2^2 = x'*(G'*G)*x
c    s.t.  sum(x)  = 1
c              x  >= 0
c
c The problem corresponds to finding the smallest vector
c (2-norm) in the convex hull of the columns of G
c
c Inputs:
c     G            -- (M x n double) matrix G, see problem above
c     varargin{1}  -- (int) maximum number of iterates allowed
c                     If not present, maxit = 100 is used
c     varargin{2}  -- (n x 1 double) vector x0 with initial (FEASIBLE) iterate.
c                     If not present, (or requirements on x0 not met) a
c                     useable default x0 will be used
c
c Outputs:
c     x       -- Optimal point attaining optimal value
c     d = G*x -- Smallest vector in the convex hull
c     q       -- Optimal value found = d'*d
c     info    -- Run data:
c                info(1) =
c                   0 = everything went well, q is optimal
c                   1 = maxit reached and final x is feasible. so q
c                       might not be optimal, but it is better than q(x0)
c                   2 = something went wrong
c                info(2) = #iterations used

subroutine qpspecial(m, n, G, maxit)
implicit none
integer          m, n, maxit
double precision G(m, n)


c  This is a program that finds the solution to the QP problem
c  []
c [x,d,q,info] = qpspecial(G,varargin)
c
c Solves the QP
c
c    min     q(x)  = || G*x ||_2^2 = x'*(G'*G)*x
c    s.t.  sum(x)  = 1
c              x  >= 0
c
c The problem corresponds to finding the smallest vector
c (2-norm) in the convex hull of the columns of G
c
c Inputs:
c     G            -- (M x n double) matrix G, see problem above
c     varargin{1}  -- (int) maximum number of iterates allowed
c                     If not present, maxit = 100 is used
c     varargin{2}  -- (n x 1 double) vector x0 with initial (FEASIBLE) iterate.
c                     If not present, (or requirements on x0 not met) a
c                     useable default x0 will be used
c
c Outputs:
c     x       -- Optimal point attaining optimal value
c     d = G*x -- Smallest vector in the convex hull
c     q       -- Optimal value found = d'*d
c     info    -- Run data:
c                info(1) =
c                   0 = everything went well, q is optimal
c                   1 = maxit reached and final x is feasible. so q
c                       might not be optimal, but it is better than q(x0)
c                   2 = something went wrong
c                info(2) = #iterations used

integer echo, idx(n)
double precision e(n), x(n)

echo = 0

c Check the dimensions

if (m*n .le. 0) then
      task = 'STOP:  Error in the dimensions of G'
      write(*,*) 'qpspecial is empty'
endif

c Create a vector of ones and use it as a starting point

do 1000 i = 1, n
      e(i) = 1d0
1000 continue

x = e

do 2100 i = 1, n
      idx(i) = (i - 1) * (n + 1) + 1
 2100  continue
