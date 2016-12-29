Program CrankNicolsonSnowline
  ! Mohamad Ali-Dib 19/01/2015
  ! Crank-Nicolson method to solve the volatiles diffusion across snowline equation. Fully working.

! Variables declaration, not all of them are used
double precision, DIMENSION(:), ALLOCATABLE :: x,s,s2,densityil2,dd,aaa,bbb,ccc,sol
Real*16, DIMENSION(:,:), ALLOCATABLE :: c
double precision :: D, deltat,deltax,deltax2,mm,nmol,mass,radius,slposition,visc,dummy,&
densityil3,nuread(200),nuread2,dummy1,tread(200),a1,a2,b1,b2,lambda,densityil,slpositionau,radiusau(200)
integer :: m,n,k,l,i,j,z,slpositionInteger,timelimit,particleind,ii,&
th2o,tco,tco2,tch4,tn2,nuind,kk, aa,bb,cc,timecounter,elapsedyears
nuread2=0

! Initialising the timecounter and Lambda (0.5 for crank-nic, 1 for implicit, 0 for explicit scheme).
timecounter=0
lambda=0.5


! Asking user to input chemical specie and the snowline position
print*, 'Please enter the particle type: 1 for water, 2 for CO, 3 for N2'
read(*,*) particleind

! This entire section is to calculate to read and calculate the average disk viscosity interior to the snowline

! This part will locate the snowline
open(unit=3, file='input.dat', action ='read')
read (3,*)
do i=1,200
read (3,*)dummy1,dummy1,dummy1,dummy1,tread(i)
if (tread(i).ge.155) then
th2o=i
elseif (tread(i).ge.50) then
tco2=i
elseif (tread(i).ge.30) then
tch4=i
elseif (tread(i).ge.25) then
tco=i
elseif (tread(i).ge.24) then
tn2=i
endif
enddo
close(3)

! This part will average the viscosity inside the snowline
open(unit=2, file='input.dat', action ='read')
read (2,*)
if (particleind.eq.1) then
do i=1,th2o
read (2,*),radiusau(i),dummy1,dummy1,nuread(i)
nuread2=nuread2+nuread(i)
enddo
nuread2=nuread2/th2o
slposition=radiusau(th2o)
elseif (particleind.eq.2) then
do i=1,tco
read (2,*),radiusau(i),dummy1,dummy1,nuread(i)
nuread2=nuread2+nuread(i)
enddo
nuread2=nuread2/tco
slposition=radiusau(tco)
elseif (particleind.eq.3) then
do i=1,tn2
read (2,*),radiusau(i),dummy1,dummy1,nuread(i)
nuread2=nuread2+nuread(i)
enddo
nuread2=nuread2/tn2
slposition=radiusau(tn2)
endif
close(2)


! Maximal time indicator allowed due to memory reasons.
timelimit=401000


! Delta_x and Delta_T
deltax=3e12
deltat=1e8


! Calculating the needed spatial grid size as a function of Delta_x and the snowline position
slpositionau=slposition/(deltax/(1.5e13))
slpositionInteger=nint(slpositionau)

! Allocating the variables sizes knowing the grid size.
Allocate(c(timelimit,slpositionInteger))
Allocate(x(slpositionInteger))
Allocate(s(slpositionInteger))
Allocate(s2(slpositionInteger))
Allocate(dd(slpositionInteger))
Allocate(aaa(slpositionInteger))
Allocate(bbb(slpositionInteger))
Allocate(ccc(slpositionInteger))
Allocate(sol(slpositionInteger))

! initilisations and initial conditions
! Source = 0
s=0

! Snowline ice density = 0
densityil=0

! Boundary values: Concentration is 1 everywhere except on the snowline
c(1,:)=1
c(1,slpositionInteger)=0


! Print spatial grid size and the snowline position
print*, slpositioninteger,slposition



open(unit=2, file='output.txt', action ='write')

! Initializing the iterating variables
i=1
j=1

! aaa,bbb and ccc are respectively the matrix lower, main and upper diagonal arrays 
! dd is the right hand term in the equation M(aa,..)X=dd.
! Boundary conditions:
bbb(1)=1
ccc(1)=-1
dd(1)=0
aaa(slpositionInteger)=0
bbb(slpositionInteger)=1
dd(slpositionInteger)=0

! ------- Main time loop -------
do while (i.le.timelimit)

! Calculating the ice density beyond the snowline by substracting the vapor concentration left from its initial value.
densityil=0
do jj=1,slpositionInteger-1
densityil=densityil+(1-c(i,jj))
enddo

! ----- Main spatial loop ------
do j=2,(slpositionInteger-1)

! Calculating the diffusivity (3*Nu to be consistent with the value used in the particle's drift module)
D=3*nuread2

! Calculating the distance 
x(j)=deltax*j

! Source function: Here comes the results of the other modules (the timescales)
if ((particleind.eq.1)) then
s((slpositionInteger-5))=(densityil*deltat)/(2e4*365*24*3600)
elseif ((particleind.ne.1)) then
s((slpositionInteger-5))=(densityil*deltat)/(3e5*365*24*3600)
endif

! Calculating the Matrix coefficients
a1=D*lambda*deltat/(deltax*x(j))
a2=D*(1-lambda)*deltat/(deltax*x(j))
b1=D*lambda*deltat/(deltax**2)
b2=D*(1-lambda)*deltat/(deltax**2)
ccc(j)=-a1-b1
bbb(j)=1+2*b1
aaa(j)=a1-b1
dd(j)=c(i,j+1)*(a2+b2)+c(i,j)*(1-2*b2)+c(i,j-1)*(-a2+b2)+s(j)

! Whatever
elapsedyears=i

! Writing the output
if ((elapsedyears.eq.4).or.(elapsedyears.eq.630).or.(elapsedyears.eq.1576)&
.or.(elapsedyears.eq.3547).or.(elapsedyears.eq.20000).or.(elapsedyears.eq.63000)&
.or.(elapsedyears.eq.250000).or.(elapsedyears.eq.390000).or.&
(elapsedyears.eq.800000).or.(elapsedyears.eq.3200000).or.(elapsedyears.eq.6400000)&
.or.(elapsedyears.eq.12800000).or.(elapsedyears.eq.42800000).or.(elapsedyears.eq.82800000)) then
if (c(i,2).ge.(1e-5)) then
if (j.eq.2) then
write(2,11), (elapsedyears*deltat)/(365*3600*24),deltax*6.6e-14,c(i,1),densityil,D
endif
write(2,11), (elapsedyears*deltat)/(365*3600*24),x(j)*6.6e-14,c(i,j),densityil,D
if (j.eq.slpositionInteger-1) then
write(2,11), (elapsedyears*deltat)/(365*3600*24),deltax*slpositionInteger*6.6e-14,c(i,slpositionInteger),densityil,D
write(2,*)
endif
endif
endif

enddo

! Calling the solver
call solve_tridiag(aaa,bbb,ccc,dd,sol,slpositionInteger)

! Reseting the problem for next iteration
do jj=1,slpositionInteger
c(i+1,jj)=sol(jj)
enddo

! Iterating
i=i+1

enddo


call system ('gnuplot gnuplotinput.plt -persist')


11 format(2x,E13.7,2x,E13.7,2x,E13.7,2x,E13.7,2x,E13.7,2x,E13.7,2x,E13.7,2x,E13.7)


End Program


Subroutine solve_tridiag(a,b,c,d,x,n)
! Thomas algorithm (simplified form of Gaussian elimination
! that can be used to solve tridiagonal systems of equations)

           implicit none
!        a - sub-diagonal (means it is the diagonal below the main diagonal)
!        b - the main diagonal
!        c - sup-diagonal (means it is the diagonal above the main diagonal)
!        d - right part
!        x - the answer
!        n - number of equations

        integer,intent(in) :: n
        real(8),dimension(n),intent(in) :: a,b,c,d
        real(8),dimension(n),intent(out) :: x
        real(8),dimension(n) :: cp,dp
        real(8) :: m
        integer i

! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i))/m
         enddo
! initialize x
         x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do
!print*, x
return
end subroutine solve_tridiag
