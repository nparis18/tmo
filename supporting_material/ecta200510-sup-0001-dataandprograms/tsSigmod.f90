module tsSigmod
	use kinkedSigmod
	implicit none
	
	real		:: global_n
	contains
	
	function gettsSig(c,s) result(val)
		real	:: c,s(:,:)
		real, allocatable	:: val(:,:)
		real	:: rho
		integer	:: i,j,n
		
		n=size(s,2)
		allocate(val(n,n))
		
		rho=exp(-c/n)
		do i=1,n
			do j=1,i
				val(i,j)=rho**abs(i-j)
				val(j,i)=val(i,j)
			enddo
		enddo
	end function	
	
	
	function gettsdelcovs(n,lam0,c0) result(val)
		use qdag_int
		real	:: lam0,c0,c1,val(0:n-1)
		integer	:: n
		real, external		:: ext_tsdelcorr
		integer	:: i
		
		global_n=n
		global_c02=c0**2
		do i=0,n-1
			global_r=i
			call qdag(ext_tsdelcorr,0.0,lam0,val(i),errabs=0.000001,irule=6)
		enddo
	end function
	
	function gettsdelSig(lam0,c0,s) result(val)
		real	:: lam0,c0,s(:,:)
		real	:: val(size(s,2),size(s,2))
		integer	:: n,i,j
		n=size(s,2)
		val=getcSig(c0,s)-toeplitz(gettsdelcovs(n,lam0,c0))
	end function

	function tsspecdens(lam,c2) result(val)
		real	:: lam,c2,val
		real	:: rho
		rho=exp(-sqrt(c2)/global_n)
		val=((1-rho**2)/Pi)/(1-2*rho*cos(lam)+rho**2)
	end function

end module

function ext_tsdelcorr(lam) result(val)
	use tsSigmod
	implicit none
	real	:: lam,val
	val=cos(lam*global_r)*tsspecdens(lam,global_c02)
end function
		