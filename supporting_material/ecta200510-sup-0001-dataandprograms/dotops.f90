module dotops
!DIR$ NOOPTIMIZE

	implicit none
	
	interface operator(.cv.)
		module procedure col_scalarr
	end interface

	interface operator(.cvc.)
		module procedure col_vecc
	end interface

	interface operator(.cvr.)
		module procedure col_vecr
	end interface

	interface operator(.clr.)
		module procedure col_leftright_v,col_leftright_vx,col_leftright_m
	end interface

	interface operator(.cud.)
		module procedure col_updown_v,col_updown_vx,col_updown_m
	end interface
	
	interface operator(.dpr.)
		module procedure dotplusr, dotplusrx
	end interface	

	interface operator(.dpc.)
		module procedure dotplusc, dotpluscx
	end interface	

	interface operator(.dmr.)
		module procedure dotminusr, dotminusrx
	end interface	

	interface operator(.dmc.)
		module procedure dotminusc, dotminuscx
	end interface

	interface operator(.dtr.)
		module procedure dottimesr, dottimesrx
	end interface	

	interface operator(.dtc.)
		module procedure dottimesc, dottimescx
	end interface	

	interface operator(.ddr.)
		module procedure dotdivider, dotdividerx
	end interface	

	interface operator(.ddc.)
		module procedure dotdividec, dotdividecx
	end interface

contains 
	function col_scalarr(a,b) result(v)
		implicit none
		real,intent(in) :: a,b
		real:: v(2)
		v(1)=a
		v(2)=b
	end function

	function col_vecr(a,b)
		implicit none
		real,intent(in) :: a(:),b(:)
		real::col_vecr(2,size(b))
		if(size(a).ne.size(b)) then
			write(*,*) 'error in col_vecr'
			stop
		endif
		
		col_vecr(1,:)=a
		col_vecr(2,:)=b
	end function

	function col_vecc(a,b)
		implicit none
		real,intent(in) :: a(:),b(:)
		real::col_vecc(size(a),2)

		if(size(a).ne.size(b)) then
			write(*,*) 'error in col_vecr'
			stop
		endif

		col_vecc(:,1)=a
		col_vecc(:,2)=b
	end function

	function col_leftright_v(A,b)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::col_leftright_v(size(A,1),size(A,2)+1)
		
		if(size(A,1).ne.size(b)) then
			write(*,*) 'error in col_leftright_v'
			stop
		endif
		col_leftright_v(:,:size(A,2))=A
		col_leftright_v(:,size(A,2)+1)=b
	end function

	function col_leftright_vx(b,A)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::col_leftright_vx(size(A,1),size(A,2)+1)
		
		if(size(A,1).ne.size(b)) then
			write(*,*) 'error in col_leftright_v'
			stop
		endif
		col_leftright_vx(:,2:)=A
		col_leftright_vx(:,1)=b
	end function
	
	function col_leftright_m(A,B)
		implicit none
		real,intent(in) :: A(:,:),B(:,:)
		real::col_leftright_m(size(A,1),size(A,2)+size(B,2))
		
		if(size(A,1).ne.size(B,1)) then
			write(*,*) 'error in col_leftright_m'
			stop
		endif
		col_leftright_m(:,:size(A,2))=A
		col_leftright_m(:,size(A,2)+1:)=B
	end function
	
	function col_updown_v(A,b)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::col_updown_v(size(A,1)+1,size(A,2))
		
		if(size(A,2).ne.size(b)) then
			write(*,*) 'error in col_updown_v'
			stop
		endif
		col_updown_v(:size(A,1),:)=A
		col_updown_v(size(A,1)+1,:)=b
	end function

	function col_updown_vx(b,A)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::col_updown_vx(size(A,1)+1,size(A,2))
		
		if(size(A,2).ne.size(b)) then
			write(*,*) 'error in col_updown_v'
			stop
		endif
		col_updown_vx(2:size(A,1)+1,:)=A
		col_updown_vx(1,:)=b
	end function
	
	function col_updown_m(A,B)
		implicit none
		real,intent(in) :: A(:,:),B(:,:)
		real::col_updown_m(size(A,1)+size(B,1),size(A,2))
		
		if(size(A,2).ne.size(B,2)) then
			write(*,*) 'error in col_updown_m'
			stop
		endif
		col_updown_m(:size(A,1),:)=A
		col_updown_m(size(A,1)+1:,:)=B
	end function
	

	function dotplusr(A,b)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::dotplusr(size(A,1),size(A,2))
		integer:: i

		if(size(A,2).ne.size(b)) then
			write(*,*) 'error in dotplusr'
			stop
		endif
		do i=1,size(A,2)
			dotplusr(:,i)=A(:,i)+b(i)
		enddo
	end function

	function dotplusrx(b,A)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::dotplusrx(size(A,1),size(A,2))
		integer:: i

		if(size(A,2).ne.size(b)) then
			write(*,*) 'error in dotplusrx'
			stop
		endif
		do i=1,size(A,2)
			dotplusrx(:,i)=A(:,i)+b(i)
		enddo
	end function

	function dotplusc(A,b)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::dotplusc(size(A,1),size(A,2))
		integer:: i

		if(size(A,1).ne.size(b)) then
			write(*,*) 'error in dotplusc'
			stop
		endif
		do i=1,size(A,2)
			dotplusc(:,i)=A(:,i)+b
		enddo
	end function

	function dotpluscx(b,A)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::dotpluscx(size(A,1),size(A,2))
		integer:: i

		if(size(A,1).ne.size(b)) then
			write(*,*) 'error in dotpluscx'
			stop
		endif
		do i=1,size(A,2)
			dotpluscx(:,i)=A(:,i)+b
		enddo
	end function
	
	function dotminusr(A,b)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::dotminusr(size(A,1),size(A,2))
		integer:: i

		if(size(A,2).ne.size(b)) then
			write(*,*) 'error in dotminusr'
			stop
		endif
		do i=1,size(A,2)
			dotminusr(:,i)=A(:,i)-b(i)
		enddo
	end function

	function dotminusrx(b,A)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::dotminusrx(size(A,1),size(A,2))
		integer:: i

		if(size(A,2).ne.size(b)) then
			write(*,*) 'error in dotminusrx'
			stop
		endif
		do i=1,size(A,2)
			dotminusrx(:,i)=b(i)-A(:,i)
		enddo
	end function

	function dotminusc(A,b)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::dotminusc(size(A,1),size(A,2))
		integer:: i

		if(size(A,1).ne.size(b)) then
			write(*,*) 'error in dotminusc'
			stop
		endif
		do i=1,size(A,2)
			dotminusc(:,i)=A(:,i)-b
		enddo
	end function

	function dotminuscx(b,A)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::dotminuscx(size(A,1),size(A,2))
		integer:: i

		if(size(A,1).ne.size(b)) then
			write(*,*) 'error in dotminuscx'
			stop
		endif
		do i=1,size(A,2)
			dotminuscx(:,i)=b-A(:,i)
		enddo
	end function

	function dottimesr(A,b)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::dottimesr(size(A,1),size(A,2))
		integer:: i

		if(size(A,2).ne.size(b)) then
			write(*,*) 'error in dottimesr'
			stop
		endif
		do i=1,size(A,2)
			dottimesr(:,i)=A(:,i)*b(i)
		enddo
	end function

	function dottimesrx(b,A)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::dottimesrx(size(A,1),size(A,2))
		integer:: i

		if(size(A,2).ne.size(b)) then
			write(*,*) 'error in dottimesrx'
			stop
		endif
		do i=1,size(A,2)
			dottimesrx(:,i)=A(:,i)*b(i)
		enddo
	end function

	function dottimesc(A,b)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::dottimesc(size(A,1),size(A,2))
		integer:: i

		if(size(A,1).ne.size(b)) then
			write(*,*) 'error in dottimesc'
			stop
		endif
		do i=1,size(A,2)
			dottimesc(:,i)=A(:,i)*b
		enddo
	end function

	function dottimescx(b,A)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::dottimescx(size(A,1),size(A,2))
		integer:: i

		if(size(A,1).ne.size(b)) then
			write(*,*) 'error in dottimescx'
			stop
		endif
		do i=1,size(A,2)
			dottimescx(:,i)=A(:,i)*b
		enddo
	end function
	
	function dotdivider(A,b)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::dotdivider(size(A,1),size(A,2))
		integer:: i

		if(size(A,2).ne.size(b)) then
			write(*,*) 'error in dotdivider'
			stop
		endif
		do i=1,size(A,2)
			dotdivider(:,i)=A(:,i)/b(i)
		enddo
	end function

	function dotdividerx(b,A)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::dotdividerx(size(A,1),size(A,2))
		integer:: i

		if(size(A,2).ne.size(b)) then
			write(*,*) 'error in dotdividerx'
			stop
		endif
		do i=1,size(A,2)
			dotdividerx(:,i)=b(i)/A(:,i)
		enddo
	end function

	function dotdividec(A,b)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::dotdividec(size(A,1),size(A,2))
		integer:: i

		if(size(A,1).ne.size(b)) then
			write(*,*) 'error in dotdividec'
			stop
		endif
		do i=1,size(A,2)
			dotdividec(:,i)=A(:,i)/b
		enddo
	end function

	function dotdividecx(b,A)
		implicit none
		real,intent(in) :: A(:,:),b(:)
		real::dotdividecx(size(A,1),size(A,2))
		integer:: i

		if(size(A,1).ne.size(b)) then
			write(*,*) 'error in dotdividecx'
			stop
		endif
		do i=1,size(A,2)
			dotdividecx(:,i)=b/A(:,i)
		enddo
	end function

end module