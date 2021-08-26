program trajectory

    implicit none

	real :: x(2400001), y(2400001), u(2400001), v(2400001)
	real :: xmid,ymid,axmid,aymid,umid,vmid
	real :: ax, ay
	real :: ms, mj, r, rj, rs, G, w
	real :: dt
	integer :: i
	
	G = 6.67e-11
	w = 1.672e-8
	ms = 1.989e30
	mj = 1.898e27
	r = 7.8e11

    dt = 3600*24*30
	
	rs = 743604142
	rj = 7.7926e11

	x(1) = rj-r/2
	y(1) = (3**0.5)*r/2
	u(1) = 0
	v(1) = 0

    do i = 1,2400000

        axmid = ax(x(i),y(i),u(i),v(i))
        aymid = ay(x(i),y(i),u(i),v(i))

        umid = u(i) + 0.5*axmid*dt
        vmid = v(i) + 0.5*aymid*dt

        xmid = x(i) + 0.5*u(i)*dt
        ymid = y(i) + 0.5*v(i)*dt

        axmid = ax(xmid,ymid,umid,vmid)
        aymid = ay(xmid,ymid,umid,vmid)

        u(i+1) = u(i) + axmid*dt
        v(i+1) = v(i) + aymid*dt

        x(i+1) = x(i) + umid*dt
        y(i+1) = y(i) + vmid*dt


    end do

    open(unit = 10, file = 'Curry.txt')

        do i = 1,2400000
            write (10,*) x(i),y(i)
        end do

    close(unit = 10)

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
