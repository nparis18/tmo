module locmod
	use globals
	implicit none
	
	integer, parameter	:: nstates=48
	
	contains
	
	function getstatelocs(istate,lightpower,n,gdens) result(val)
		integer	:: istate,n
		real	:: lightpower,val(2,n)
		real, optional	:: gdens(n)
		
		logical, save	:: initflag=.true.		
		real, allocatable, save	:: mdata(:,:)
		character(len=:),allocatable,save	:: names(:)		
		real, allocatable	:: sdata(:,:)
		real, allocatable	:: dG(:),G(:)
		
		real	:: u(3),s0(2),s(2)
		integer	:: i,j,l
		
		if(initflag) then
			mdata=loadcsv("c:/dropbox/SHAR/data/US_Light_ext.csv")
			names=loadstrings("c:/dropbox/SHAR/data/US_Light_extnames.txt")
			initflag=.false.
		endif
		if(istate>nstates) then
			print *,"getstatelocs encountered state with no data"
			stop
		endif
!		print *,"generating random data from state",istate, names(istate:istate)
		sdata=selectifr(mdata,mdata(:,1)==istate)

		dG=sdata(:,6)
		dG=merge(dG,0.0,isfinite(dG))+1E-100
		dG=dG**lightpower
		G=dG
		do i=2,size(G)
			G(i)=G(i-1)+G(i)
		enddo
		G=G/(0.99999999*G(size(G)))
		do i=1,n
			call rnun(u)
			j=count(u(3)>G)+1		
			val(:,i)=sdata(j,2:3)+u(1:2)*(sdata(j,4:5)-sdata(j,2:3))
			if(present(gdens)) gdens(i)=dG(j)
		enddo
	end function
			

end module