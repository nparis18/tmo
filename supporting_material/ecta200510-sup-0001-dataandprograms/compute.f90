module compute	
	use locmod
	use rpmod
	use rootfinder
	use Sigmod
	implicit none
	
	contains
	

	function getcvfromOms(Oms,q,level) result(val)
		real	:: Oms(:,:,:),level,val
		integer	:: q
		real	:: cv0,cv1,cv,rps(size(Oms,3))
		integer	:: i,j
		
		i=1
		cv0=getStudentcv(q,level)
		do
			cv1=cv0
			do
				if(getrp(Oms(1:q+1,1:q+1,i),cv1)>level) then
					cv0=cv1
					cv1=cv1+1.0/sqrt(real(q))
				else
					exit
				endif
			enddo
			do 
				if(cv1-cv0<0.001/sqrt(real(q))) exit
				cv=(cv0+cv1)/2
				if(getrp(Oms(1:q+1,1:q+1,i),cv)>level) then
					cv0=cv
				else
					cv1=cv
				endif
			enddo
			do j=1,size(rps)
				rps(j)=getrp(Oms(1:q+1,1:q+1,j),cv1)
			enddo
			i=maxloc(rps,dim=1)
			cv0=cv1
			if(rps(i)<level) exit
		enddo
		val=cv1
	end function
	
	function getOms(W,s,c0) result(Oms)
		real	:: W(:,:),s(:,:),c0
		real	:: Oms(size(W,2),size(W,2),nc)
		integer				:: i
		real, parameter	:: cfacs(nc)=exp([(5.0*(i-1.0)/(nc-1),i=1,nc)])
		real, allocatable	:: Sig(:,:),Om(:,:)
		
		do i=1,nc
			Sig=getcSig(c0*cfacs(i),s)
			Oms(:,:,i)=matmul(transpose(W),matmul(Sig,W))
		enddo
	end function
			
	subroutine setSCPCqcv(W,s,c0,q,cv)
		real, allocatable	:: W(:,:)
		real	:: s(:,:),c0,cv
		integer	:: q
		real	:: Oms(size(W,2),size(W,2),nc)
		real	:: cvs(size(W,2)-1),lengths(size(W,2)-1)

		Oms=getOms(W,s,c0)
		do q=1,size(W,2)-1
			cvs(q)=getcvfromOms(Oms,q,level)	
			lengths(q)=getIlength(q,cvs(q))
		enddo
		q=minloc(lengths,1)
		cv=cvs(q)
		W=W(:,1:1+q)
	end subroutine
		
	subroutine setSCPCqcv_alt(W,s,Om1,c0,q,cv)
		real, allocatable	:: W(:,:)
		real	:: s(:,:),c0,cv,Om1(:,:)
		integer	:: q
		real	:: Oms(size(W,2),size(W,2),nc)
		real	:: cvs(size(W,2)-1),lengths(size(W,2)-1)
		Oms=getOms(W,s,c0)
		do q=1,size(W,2)-1
			cvs(q)=getcvfromOms(Oms,q,level)	
			lengths(q)=getlength(Om1(1:q+1,1:q+1),cvs(q))
		enddo
		q=minloc(lengths,1)
		cv=cvs(q)
		W=W(:,1:1+q)
	end subroutine
	
end module


	