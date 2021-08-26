program trajectory

    implicit none

	real :: x(520001), y(520001), u(520001), v(520001), p(520000)
	real :: dxdt,dydt,dudt,dvdt
	real :: kx1,ky1,ku1,kv1,kx2,ky2,ku2,kv2,kx3,ky3,ku3,kv3,kx4,ky4,ku4,kv4
	real :: ms, mj, r, rj, rs, w, G
	real :: dt
	integer :: i
	
	G = 6.67e-11
	ms = 1.989e30
	mj = 1.898e27
	r = 7.8e11

    dt = 3600*24*7
	
	rs = 743604142
	rj = 7.7926e11
	
	w = 1.68e-8

	x(1) = rj - r/2
	y(1) = (3**0.5)*r/2
	u(1) = 0
	v(1) = 0

    do i = 1,520000

        kx1 = dxdt(x(i),y(i),u(i),v(i))
		ky1 = dydt(x(i),y(i),u(i),v(i))
		ku1 = dudt(x(i),y(i),u(i),v(i))
		kv1 = dvdt(x(i),y(i),u(i),v(i))
		
		kx2 = dxdt(x(i)+dt*kx1/2,y(i)+dt*ky1/2,u(i)+dt*ku1/2,v(i)+dt*kv1/2)
		ky2 = dydt(x(i)+dt*kx1/2,y(i)+dt*ky1/2,u(i)+dt*ku1/2,v(i)+dt*kv1/2)
		ku2 = dudt(x(i)+dt*kx1/2,y(i)+dt*ky1/2,u(i)+dt*ku1/2,v(i)+dt*kv1/2)
		kv2 = dvdt(x(i)+dt*kx1/2,y(i)+dt*ky1/2,u(i)+dt*ku1/2,v(i)+dt*kv1/2)
		
		kx3 = dxdt(x(i)+dt*kx2/2,y(i)+dt*ky2/2,u(i)+dt*ku2/2,v(i)+dt*kv2/2)
		ky3 = dydt(x(i)+dt*kx2/2,y(i)+dt*ky2/2,u(i)+dt*ku2/2,v(i)+dt*kv2/2)
		ku3 = dudt(x(i)+dt*kx2/2,y(i)+dt*ky2/2,u(i)+dt*ku2/2,v(i)+dt*kv2/2)
		kv3 = dvdt(x(i)+dt*kx2/2,y(i)+dt*ky2/2,u(i)+dt*ku2/2,v(i)+dt*kv2/2)
		
		kx4 = dxdt(x(i)+dt*kx3,y(i)+dt*ky3,u(i)+dt*ku3,v(i)+dt*kv3)
		ky4 = dydt(x(i)+dt*kx3,y(i)+dt*ky3,u(i)+dt*ku3,v(i)+dt*kv3)
		ku4 = dudt(x(i)+dt*kx3,y(i)+dt*ky3,u(i)+dt*ku3,v(i)+dt*kv3)
		kv4 = dvdt(x(i)+dt*kx3,y(i)+dt*ky3,u(i)+dt*ku3,v(i)+dt*kv3)
		
		x(i+1) = x(i) + 0.16667*dt*(kx1+2*kx2+2*kx3+kx4)
		y(i+1) = y(i) + 0.16667*dt*(ky1+2*ky2+2*ky3+ky4)
		u(i+1) = u(i) + 0.16667*dt*(ku1+2*ku2+2*ku3+ku4)
		v(i+1) = v(i) + 0.16667*dt*(kv1+2*kv2+2*kv3+kv4)
		
		p(i) = -0.5*(w**2)*(x(i+1)**2+y(i+1)**2) - (G*ms)/(sqrt((x(i+1)+rs)**2+y(i+1)**2)) - (G*mj)/(sqrt((x(i+1)-rj)**2+y(i+1)**2))

    end do

    open(unit = 10, file = 'RKP.txt')

        do i = 1,520000
            write (10,*) i, p(i)
        end do

    close(unit = 10)
	
	print *, 2*(maxval(p)-minval(p))/(maxval(p)+minval(p))

end program

real function dxdt(x,y,u,v)

	real :: x,y,u,v

	dxdt = u

end function

real function dydt(x,y,u,v)

	real :: x,y,u,v

	dxdt = v

end function

real function dudt(x,y,u,v)

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


	dudt = 2*w*v + (w**2)*x - (G*ms*(x+rs))/(((x+rs)**2+y**2)**1.5) - (G*mj*(x-rj))/(((x-rj)**2+y**2)**1.5)

end function

real function dvdt(x,y,u,v)

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


	dvdt = -2*w*u + (w**2)*y - (G*ms*y)/(((x+rs)**2+y**2)**1.5) - (G*mj*y)/(((x-rj)**2+y**2)**1.5)

end function
