program potential 

	implicit none
	
	real, dimension (2000) :: x, y
	integer :: i, j, num
	real :: n
	real :: pot
	real :: lim
	
	lim = 0.25
	n = size(x)
	num = 2000
	
	do i = 1,num

		x(i) = -2 + (4/n)*i
		y(i) = -2 + (4/n)*i
	
	end do

	open(unit = 10, file='potential.txt')

		do i = 1,num
		
			do j = 1,num
			
				if (pot(x(i),y(j))<lim) then
	
					write(10,*) x(i), y(j), pot(x(i),y(j))
				
				else 
					
					write(10,*) x(i), y(j), 0
				
				end if
			 
			end do
			
		end do

end program

real function pot(x,y)
	
	real :: x,y
	real, parameter :: ms = 10, mj = 1, G = 1, w = 3.16, rs = .09, rj = .91
	
	potential =  (0.5*(w**2)*(x**2+y**2) + (G*ms/sqrt((x+rs)**2 + y**2)) + (G*mj/sqrt((x-rj)**2 + y**2)))/100
	
end function