module stepSigmod
	use kinkedSigmod
	implicit none
	
	contains
	
	function getstepcorr(r,lam0,c0,initflag) result(val)
		use qdag_int
		use qdval_int
		real	:: r,lam0,c0,c1,val
		logical				:: initflag
		integer				:: i
		integer, parameter	:: nr=100
		real, parameter		:: rgrid(nr)=[(sqrt(2.0)*((i-1.0)/(nr-1))**2,i=1,nr)]
		real, save			:: fgrid(nr)
		real, external		:: ext_stepcorr
		if(initflag) then
			global_c02=c0**2
			do i=1,nr
				global_r=rgrid(i)
				call qdag(ext_stepcorr,0.0,lam0,fgrid(i),errabs=0.000001,irule=6)
			enddo
			val=-1	
			return
		endif
		val=qdval(r,rgrid,fgrid,check=.false.)
	end function
	
	
	function getstepSig(lam0,c0,s) result(val)
		real	:: lam0,c0,c1,s(:,:)
		real, allocatable	:: val(:,:)
		integer	:: n,i,j
		n=size(s,2)
		allocate(val(n,n))
		val(1,1)=getstepcorr(0.0,lam0,c0,.true.)
!$omp parallel do private(j)
		do i=1,n
			do j=1,i
				val(j,i)=getstepcorr(norm2(s(:,i)-s(:,j)),lam0,c0,.false.)
				val(i,j)=val(j,i)
			enddo
		enddo
		val=getcSig(c0,s)-val
	end function
	
end module

function ext_stepcorr(lam) result(val)
	use kinkedSigmod
	implicit none
	real	:: lam,val
	val=lam*Bessel_J0(lam*global_r)*specdens(lam,global_c02)
end function

	
	