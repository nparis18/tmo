module robmod
	use tsSigmod
	implicit none
	
	
	contains

	function getAfromOm(Om) result(val)
		use evasf_int
		real	:: Om(:,:),val(size(Om,1),size(Om,1))
		real	:: chol(size(Om,1),size(Om,1)),cholx(size(Om,1),size(Om,1)),evals(1)
		chol=choleski(Om)
		cholx(1,:)=chol(1,:)
		cholx(2:,:)=-chol(2:,:)
		call evasf(1,matmul(transpose(chol),cholx),.false.,evals)
		val(:,1)=Om(:,1)/evals(1)
		val(:,2:)=-Om(:,2:)/evals(1)
	end function

	function getnu(A,lamA0,P,Pinv) result(val)
		use evlsf_int
		real	:: A(:,:),lamA0(:),P(:,:),Pinv(:,:),val(size(A,1)-1)
		real	:: B(size(A,1),size(A,1)),lamB(size(A,1)),nu(size(A,1)-1)
		integer	:: i

		B=matmul(Pinv,matmul(A,P))
		B=.5*(B+transpose(B))
		call evlsf(B,lamB)
		lamB=-sort(-lamB)
		nu=-lamB(2:)+lamA0(2:)*lamB(1)
		nu(1)=nu(1)+1-lamB(1)
		do i=2,size(A,1)-1
			nu(i)=nu(i-1)+nu(i)
		enddo
		val=nu
	end function
	
	subroutine initA0(Om,lamA0,P,Pinv)
		use evcrg_int
		real	:: Om(:,:),lamA0(:),P(:,:),Pinv(:,:)
		complex	:: evals(size(Om,1)),evecs(size(Om,1),size(Om,1))
		integer	:: ind(size(Om,1))

		call evcrg(getAfromOm(Om),evals,evecs)
		lamA0=-evals
		ind=sortind(lamA0)
		lamA0=-lamA0(ind)
		P=evecs(:,ind)
		Pinv=invertgen(P)
	end subroutine

	function getrobust_ts(W,c0,s) result(val)
		real	:: W(:,:),c0,s(:,:),val
		integer, parameter	:: n0=50
		real, dimension(size(W,2),size(W,2)):: Pmat,Pinv,A
		real, dimension(size(W,2))			:: lamA0
		real	:: nus(size(W,2)-1,n0)
		real	:: c1,lam0
		integer	:: i,j
		real, allocatable	:: Om1(:,:)

		lam0=Pi*min(2*c0/size(s,2),1.0)	
		c1=1.03*c0
		Om1=getOmfromW(W,gettsSig(c1,s))
		call initA0(getOmfromW(W(:,1).clr.getcv(Om1,level)*W(:,2:),gettsSig(c1,s)),lamA0,Pmat,Pinv)
		do i=1,n0
			A=getAfromOm(Om1+5*getOmfromW(W,gettsdelSig((real(i)/n0)*lam0,c1,s)))
			nus(:,i)=getnu(A,lamA0,Pmat,Pinv)
		enddo
		val=minval(nus)
	end function

	function getrobust_step(W,c0,s) result(val)
		use stepSigmod
		real	:: W(:,:),c0,s(:,:),val
		integer, parameter	:: n0=50
		real, dimension(size(W,2),size(W,2)):: Pmat,Pinv,A
		real, dimension(size(W,2))			:: lamA0
		real	:: nus(size(W,2)-1,n0),nu0(size(W,2)-1),minratio
		real	:: lam0, c1
		integer	:: i,j,failflag
		real, allocatable	:: Om1(:,:)

		lam0=7*c0
		j=1
		do
			c1=(1+.02*(j-1))*c0
			Om1=getOmfromW(W,getcSig(c1,s))
			call initA0(getOmfromW(W(:,1).clr.getcv(Om1,level)*W(:,2:),getcSig(c1,s)),lamA0,Pmat,Pinv)
			do i=1,n0
				A=getAfromOm(Om1+10*getOmfromW(W,getstepSig((real(i)/n0)*lam0,c1,s)))
				nus(:,i)=getnu(A,lamA0,Pmat,Pinv)
			enddo
			if(minval(nus)>0) exit
			j=j+1
			if(j>20) exit
		enddo
		if(j>20) then
			val=1E10
		else
			val=c1/c0
		endif
	end function
end module


	