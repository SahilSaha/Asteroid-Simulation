program trajectory

    implicit none

	real :: x(52001), y(52001), u(52001), v(52001), p(52000)
	real :: ax, ay
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

    do i = 1,52000

        x(i+1) = x(i) + u(i)*dt
		y(i+1) = y(i) + v(i)*dt
		
		u(i+1) = u(i) + ax(x(i),y(i),u(i),v(i))*dt
		v(i+1) = v(i) + ay(x(i),y(i),u(i),v(i))*dt
		
		p(i) = -0.5*(w**2)*(x(i)**2+y(i)**2) - (G*ms)/(sqrt((x(i)+rs)**2+y(i)**2)) - (G*mj)/(sqrt((x(i)-rj)**2+y(i)**2))

    end do

    open(unit = 10, file = 'EulerP.txt')

        do i = 1,52000
            write (10,*) i, p(i)
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
