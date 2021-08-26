program trajectory
	
	use omp_lib
	
    implicit none

	real, dimension(100000) :: x,y,u,v,xmid,ymid,axmid,aymid,umid,vmid,s
	real :: ax, ay
	real :: dt
	integer :: i,j
	
    dt = 3600*24*14
	
	open(unit = 20, file = 'GeneratePos.txt')
	open(unit = 30, file = 'GenerateVel.txt')
	
	do i = 1,100000
	
		read(20,*) x(i),y(i)
		read(30,*) u(i),v(i)
		s(i) = 0
	
	end do
	
	close(unit = 20)
	close(unit = 30)
	
	
	open(unit = 40, file = 'DistSmall.txt')
	open(unit = 10, file = 'AstPositSmall.txt')
	open(unit = 50, file = 'avg.txt')

	!$OMP PARALLEL DO
    do i = 1,100000
	
		do j = 1, 14000
			
			if (sqrt(x(i)**2+y(i)**2)<8e11) then

				axmid(i) = ax(x(i),y(i),u(i),v(i))
				aymid(i) = ay(x(i),y(i),u(i),v(i))

				umid(i) = u(i) + 0.5*axmid(i)*dt
				vmid(i) = v(i) + 0.5*aymid(i)*dt

				xmid(i) = x(i) + 0.5*u(i)*dt
				ymid(i) = y(i) + 0.5*v(i)*dt

				axmid(i) = ax(xmid(i),ymid(i),umid(i),vmid(i))
				aymid(i) = ay(xmid(i),ymid(i),umid(i),vmid(i))

				u(i) = u(i) + axmid(i)*dt
				v(i) = v(i) + aymid(i)*dt

				x(i) = x(i) + umid(i)*dt
				y(i) = y(i) + vmid(i)*dt
				
				s(i) = s(i) + sqrt(x(i)**2+y(i)**2)
					
			else
				
				exit
					
			end if

		end do
			
		if (sqrt(x(i)**2+y(i)**2)<8e11) then
				
			write(40,*) sqrt(x(i)**2+y(i)**2)
			write(10,*) x(i), y(i)
			write(50,*) s(i)/14000
			
		end if

    end do
	!$OMP END PARALLEL DO
	
	close(unit = 40)
	close(unit = 10)

	
	print *, "Done!"

end program

real function ax(x,y,u,v)

	real :: x,y,u,v
	real :: ms, mj, r, G, w
	real :: rj, rs
	
	ms = 1.989e30
	mj = 1.898e27
	r = 7.8e11
	G = 6.67e-11
	w = 1.68e-8

	rs = 743604142
	rj = 7.7926e11


	ax = 2*w*v + (w**2)*x - (G*ms*(x+rs))/(((x+rs)**2+y**2)**1.5) - (G*mj*(x-rj))/(((x-rj)**2+y**2)**1.5)

end function

real function ay(x,y,u,v)

	real :: x,y,u,v
	real :: ms , mj, r, G, w
	real :: rj, rs
	
	ms = 1.989e30
	mj = 1.898e27
	r = 7.8e11
	G = 6.67e-11
	w = 1.68e-8

	rs = 743604142
	rj = 7.7926e11


	ay = -2*w*u + (w**2)*y - (G*ms*y)/(((x+rs)**2+y**2)**1.5) - (G*mj*y)/(((x-rj)**2+y**2)**1.5)

end function
