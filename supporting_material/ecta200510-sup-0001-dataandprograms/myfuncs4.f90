!INCLUDE 'link_fnl_shared_imsl.h'    
INCLUDE 'link_fnl_static_imsl.h'    

module myfuncs
	use linds_int
	use lftds_int
	use lfdds_int

	use rnopt_int
	use rnnoa_int
	use rnchi_int
	use rnunf_int
	use rnund_int
	use rnstt_int
	use rnset_int
	use rnfdf_int
	use rnun_int
	
	use tdf_int
	use tpr_int
	use tin_int
	use fin_int
	use anorin_int
	use anordf_int
	use anorpr_int
	use eqtil_int
	use ordst_int
	use svrgn_int
	use svrgp_int

	use uminf_int
	use bconf_int
	use erset_int
	use u4inf_int
!    use neqnf_int
	
	use qdval_int
	use bsnak_int
	use bscpp_int
	use csdec_int
	use pp1gd_int
	use bsitg_int
	use bs1gd_int
	
	use gamit_int
	use bsi0_int
	use bsi0e_int

	use IERCD_int
	
	use gqrul_int
	
!	use ML
!	use dispmodule

	implicit none
!	use table_io
!	use imsl_libraries
	integer(8)	time0, time1
	real, parameter		:: pi=3.14159265358979323846
	
!	private
	
!	public rann,invertpd,choleski,detpd,logdetpd,ransubsample,logmeanexp

    interface rann
        module procedure rann_m, rann_v, rann_s
	end interface

    interface setbounds
        module procedure setbounds_m,setbounds_v
	end interface

    interface mdisp
        module procedure mydispvec,mydispmat,mydispscalar,mydispstring
	end interface

    interface quantile
        module procedure quantile_s,quantile_v
	end interface
	
	interface gausscdfinv
		module procedure gausscdfinv_s, gausscdfinv_v
	end interface
    
    contains

	subroutine mkGQxw(GQxw)
		implicit none
		
		real	:: GQxw(:,:)

		call gqrul(size(GQxw,1),GQxw(:,1),GQxw(:,2))
		GQxw(:,1)=.5+.5*GQxw(:,1)
		GQxw(:,2)=.5*GQxw(:,2)
	end subroutine
	
	function gausscdfinv_s(x) result(val)
		implicit none
		real	:: x,val,a(1),y(1)
		a=x
		call vdcdfnorminv( 1, a, y )
		val=y(1)
	end function
	
	function gausscdfinv_v(x) result(val)
		implicit none
		real	:: x(:),val(size(x))
		call vdcdfnorminv( size(x), x, val )
	end function

	function gausscdf(x) result(val)
		implicit none
		real	:: x,val,a(1),y(1)
		a=x
		call vdcdfnorm( 1, a, y )
		val=y(1)
	end function
	
	function gausscdfvec(x) result(val)
		implicit none
		real	:: x(:),val(size(x))
		call vdcdfnorm( size(x), x, val )
	end function
	
	elemental function gaussdens(x) result (val)
		implicit none
		real, intent(in)	:: x
		real				:: val
		val=(2*pi)**(-0.5)*exp(-.5*x**2)
	end function
	
	elemental function gausspdf(x) result (val)
		implicit none
		real, intent(in)	:: x
		real				:: val
		val=gaussdens(x)
	end function

	function rotmat(phi) result(mat)
		real	:: phi,mat(2,2)
		mat=reshape([cos(phi),sin(phi),-sin(phi),cos(phi)],[2,2])
	end function
		
	subroutine setbounds_m(m,l,u)
		implicit none
		real	:: m(:,:),l,u
		where(m<l) 
			m=l
		elsewhere	
			where(m>u) m=u
		endwhere
	end subroutine
	
	subroutine setbounds_v(v,l,u)
		implicit none
		real	:: v(:),l,u
		where(v<l) 
			v=l
		elsewhere	
			where(v>u) v=u
		endwhere
	end subroutine

	subroutine inittime()
		call system_clock(time0)
	end subroutine
	
	subroutine printtime()
!		character(len=1024) output
		integer(8), parameter	:: nrate=1000
		integer(8)	:: h,m,s,t
		call system_clock(time1)
		t=(time1-time0)/1000
		h=t/(3600*nrate)
		t=t-(3600*nrate)*h
		m=t/(nrate*60)
		t=t-(nrate*60)*m
		s=t/nrate
		t=t-nrate*s
	
		print *,"time elapsed ",h,":",m,":",s,".",t
!		call VSprint(output)
	end subroutine
	
	elemental function isfinite(x) result(val)
		use, intrinsic :: ieee_arithmetic 
		real, intent(in)	:: x
		logical	:: val
		val =ieee_is_finite(x)
	end function
		
	elemental function boole(flag) result(val)
		logical, intent(in)	:: flag
		integer	:: val
		if(flag) then 
			val=1
		else
			val=0
		endif
	end function
	
	elemental function rboole(flag) result(val)
		logical, intent(in)	:: flag
		real	:: val
		if(flag) then 
			val=1
		else
			val=0
		endif
	end function

	!function itostring(i) result(val)
	!	implicit none
	!	integer		:: i
	!	character(len=:),allocatable :: val
	!	character(len=100)	:: string
	!	write (string,*) i
	!	allocate(character(len=len(trim(adjustl(string))))::val)
	!	val=trim(adjustl(string))
	!end function
		
    function getquadmin(x,y) result(val)
        implicit none
        real    :: x(3),y(3),val
        real    :: dy(3),den
        
        dy(1)=y(2)-y(1)
        dy(2)=y(3)-y(1)
        dy(3)=y(3)-y(2)
        den=(x(3)*dy(1)-x(2)*dy(2)+x(1)*dy(3))
		if(den/=0) then
			val=.5*(x(3)**2*dy(1)-x(2)**2*dy(2)+x(1)**2*dy(3))/den
		else
			val=x(2)
		endif
    end function
    
    function getquadint(x,y,x0) result(val)
        implicit none
        real    :: x(3),x0,val
		real	:: y(3)
        
        val=((x0 - x(3))*((x0 - x(2))*(x(2) - x(3))*y(1) + (x0 - x(1))*(-x(1) + x(3))*y(2)) + (x0 - x(1))*(x0 - x(2))*(x(1) - x(2))*y(3))/((x(1) - x(2))*(x(1) - x(3))*(x(2) - x(3)))
        
    end function
    
    function getlinint(x,y,x0) result(val)
        implicit none
        real    :: x(2),y(2),x0,val
        
        val=((x0-x(2))*y(1)+(x(1)-x0)*y(2))/(x(1) - x(2))
        
	end function

	elemental function logit(x)
		real, intent (in)	:: x
		real				:: logit
		if(x>0) then
			logit=1.0/(1+exp(-x))
		else
			logit=exp(x)/(1+exp(x))
		endif
	end function

	elemental function invlogit(x) result(val)
		real, intent (in)	:: x
		real				:: val
		val=log(x/(1-x))
	end function

	function invertpd(A) result(val)
		real::A(:,:)
		real::val(size(A,1),size(A,2))
		call linds(A,val)
    end function
    
    function invertgen(A) result(val)
        use linrg_int
		real::A(:,:)
		real::val(size(A,1),size(A,2))
		call linrg(A,val)
	end function
	
	function conpd(A) result(val)
		use lfcds_int
		real::A(:,:),val
		real::B(size(A,1),size(A,2))
		call lfcds(A,B,val)
		val=1/val
	end function
	
	function choleski(A) result(val)
		real::A(:,:)
		real::val(size(A,1),size(A,2))
		integer	:: j1,j2
		call lftds(A,val)
		do j1=1,size(A,1)
			do j2=j1+1,size(A,1)
				val(j1,j2)=0
			enddo
		enddo
	end function
	
	function choleskird(A) result(val)
		use lchrg_int
		real::A(:,:)
		real::val(size(A,1),size(A,2))
		integer	:: j1,j2
		call lchrg(A,val)
		do j1=1,size(A,1)
			do j2=j1+1,size(A,1)
				val(j1,j2)=0
			enddo
		enddo
    end function
    
    function outerprod(v1,v2) result(val)
        real    :: v1(:),v2(:),val(size(v1),size(v2))
        integer :: i1,i2
        do i1=1,size(v1)
            do i2=1,size(v2)
                val(i1,i2)=v1(i1)*v2(i2)
            enddo
        enddo
    end function
	
	function diagonal(A) result(val)
		real::A(:,:),val(size(A,1))
		integer :: i
		do i=1,size(A,1)
			val(i)=A(i,i)
		enddo
	end function
	
	function eye(k) result(val)
		integer	:: k
		real	:: val(k,k)
		integer	:: i
		val=0
		do i=1,k
			val(i,i)=1.0
		enddo
	end function
	
	function toeplitz(v) result(val)
		real	:: v(:),val(size(v),size(v))
		integer	:: i,j
		do i=1,size(v)
			do j=1,i
				val(j,i)=v(i-j+1)
				val(i,j)=val(j,i)
			enddo
		enddo
	end function
		
	function zeros(k) result(val)
		integer	:: k
		real	:: val(k)
		val=0
	end function

	function ones(k) result(val)
		integer	:: k
		real	:: val(k)
		val=1
	end function

	function diagmat(vec) result(val)
		implicit none
		real	:: vec(:)
		real	:: val(size(vec),size(vec))
		integer	:: i
		val=0
		do i=1,size(vec)
			val(i,i)=vec(i)
		enddo
	end function
	
	function bandmat(vec) result(val)
		real	:: vec(:)
		real	:: val(size(vec)/2+1,size(vec)/2+1)
		integer	:: i
		val=0
		do i=1,size(vec)/2+1
			val(:,i)=vec(size(vec)/2+2-i:size(vec)/2+1-i+size(vec)/2)
		enddo
	end function
	
	function getdiagonal(mat) result(val)
		real	:: mat(:,:)
		real	:: val(size(mat,1))
		integer	:: i
		val=0
		do i=1,size(val)
			val(i)=mat(i,i)
		enddo
	end function

	function detpd(A)
		real::A(:,:)
		real::detpd
		real::F(size(A,1),size(A,2)),det1,det2
		
		call lftds(A,F)
		call lfdds(F,det1,det2)
		detpd=det1*10**det2
	end function

	function logdetpd(A)
		real::A(:,:)
		real::logdetpd
		real::F(size(A,1),size(A,2)),det1,det2
		
		call lftds(A,F)
		call lfdds(F,det1,det2)
		logdetpd=log(det1)+log(10.0)*det2
	end function
	
	function detgen(A) result(val)
		use lftrg_int
		use lfdrg_int
		real:: A(:,:), val
		real::F(size(A,1),size(A,2)),det1,det2
		integer	:: ipvt(size(A,1))
		
		call lftrg(A,F,ipvt)
		call lfdrg(F,ipvt,det1,det2)
		val=det1*10**det2
	end function
	
	function logabsdetgen(A) result(val)
		use lftrg_int
		use lfdrg_int
		real:: A(:,:), val
		real::F(size(A,1),size(A,2)),det1,det2
		integer	:: ipvt(size(A,1))
		
		call lftrg(A,F,ipvt)
		call lfdrg(F,ipvt,det1,det2)
		val=log(abs(det1))+log(10.0)*det2
	end function
		
	function rann_m(n,m)
		integer::n,m,i
		real::rann_m(n,m)
		
		do i=1,m
			call rnnoa(rann_m(:,i))
		enddo
	end function		
		
	function rann_v(n)
		integer::n
		real::rann_v(n)
	
		call rnnoa(rann_v)
	end function

	function rann_s()
		real::rann_s,dummy(1)
		call rnnoa(dummy)
		rann_s=dummy(1)
	end function

	function ranu_s() result(out)
		real::out
		out=rnunf()
	end function
	
	function ranu_v(n) result(out)
		integer::n
		real::out(n)
		call rnun(out)
	end function
	
	function ranui_v(n,k) result(v)
		integer		::n,k
		integer		::v(k)
		call rnund(n,v)
	end function
	
	function ranui_s(n) result(v)
		integer		::v,dummy(1)
		integer		::n
		call rnund(n,dummy)
		v=dummy(1)
	end function

	function logmeanexp(v) result(out)
		real		::v(:)
		real		::out,maxv
		maxv=maxval(v)
		out=maxv+log(sum(exp(v-maxv)))-log(real(size(v)))
	end function

	function logsumexp(v) result(out)
		real		::v(:)
		real		::out,maxv
		maxv=maxval(v)
		out=maxv+log(sum(exp(v-maxv)))
    end function
    
    function orderstat(k,a) result(val)
        integer  :: k
        real     :: a(:),val,hope,pesty
        integer  ::  n
        integer  :: l,r,l2,r2
        n=size(a)
        l = 1
        r = n
        do while (l < r)
            hope = a(k)
            l2 = l
            r2 = r
            do while (l2 <= r2)
                do while (a(l2)< hope)
                    l2 = l2 + 1
                end do
                do while (hope<a(r2))
                    r2 = r2 - 1
                end do
                if (l2 <r2) then
                    pesty = a(l2)
                    a(l2) = a(r2)
                    a(r2) = pesty
                endif
                if (l2 <= r2) then
                    l2 = l2 + 1
                    r2 = r2 - 1
                endif
            end do
            if (r2<k) l = l2
            if (k<l2) r = r2
        end do
        val = a(k)
    end function 

    function quantile_v(X,ps) result(qs)
		real		::X(:),ps(:),qs(size(ps))
        real, allocatable   :: xcopy(:)
        integer     :: i,k
        
		allocate(xcopy(size(x)+2))
		xcopy(1)=-huge(1.0); xcopy(size(x)+2)=huge(1.0)
		do i=1,size(x)
			xcopy(1+i)=x(i)
		enddo
        do i=1,size(ps)
            k=nint(ps(i)*size(x))+1
            qs(i)=orderstat(k,xcopy)
        enddo
    end function
    
    function quantile_s(X,p) result(q)
		real		::X(:),p,q
        real        :: qs(1)
        qs=quantile_v(X,[p])
        q=qs(1)
    end function

    function sort(X) result(val)
        real    :: X(:),val(size(X))
        call svrgn(x,val)
	end function
	
	function sortind(X) result(val)
		use svrgp_int
        real    :: X(:),xp(size(x))
		integer	:: val(size(X)),i
		val=[(i,i=1,size(X))]
        call svrgp(x,xp,val)
    end function
    
	function stddev(X) result(std)
		real		::X(:),std
        std=sqrt(sum((X-sum(X)/size(X))**2)/size(X))
	end function
	
	function getcorr(X,Y) result(val)
		real		::X(:),Y(:),val
		real		::xbar,ybar
		integer		:: n
		n=size(x)
		xbar=sum(x)/n
		ybar=sum(y)/n
		val=(sum(x*y)-n*xbar*ybar)/sqrt((sum(x**2)-n*xbar**2)*(sum(y**2)-n*ybar**2))
	end function
    
    function selectifc(X,cond) result(val)
		implicit none
		real	:: X(:,:)
		logical	:: cond(:)
		integer	:: i,c
		real, allocatable	:: val(:,:)
		if(size(X,2).ne.size(cond)) then
			print *,"vector size doesn't match in selectifc"
			stop
		endif
		allocate(val(size(X,1),count(cond)))
		c=1
		do i=1,size(X,2)
			if(cond(i)) then 
				val(:,c)=X(:,i)
				c=c+1
			endif
		enddo
	end function

	function selectifr(X,cond) result(val)
		implicit none
		real	:: X(:,:)
		logical	:: cond(:)
		integer	:: i,c
		real, allocatable	:: val(:,:)
		if(size(X,1).ne.size(cond)) then
			print *,"vector size doesn't match in selectifr"
			stop
		endif
		allocate(val(count(cond),size(X,2)))
		c=1
		do i=1,size(X,1)
			if(cond(i)) then 
				val(c,:)=X(i,:)
				c=c+1
			endif
		enddo
	end function
	
	function selectif_r(v,cond) result(val)
		implicit none
		real	:: v(:)
		logical	:: cond(:)
		integer	:: i,c
		real, allocatable	:: val(:)
		if(size(v).ne.size(cond)) then
			print *,"vector size doesn't match in selectif"
			stop
		endif
		allocate(val(count(cond)))
		c=1
		do i=1,size(v)
			if(cond(i)) then 
				val(c)=v(i)
				c=c+1
			endif
		enddo
	end function
	
	function reversec(A) result(val)
		implicit none
		real	:: A(:,:)
		real	:: val(size(A,1),size(A,2))
		integer	:: i,j
		val=A([(j,j=size(A,1),1,-1)],:)
	end function
	
	function reverser(A) result(val)
		implicit none
		real	:: A(:,:)
		real	:: val(size(A,1),size(A,2))
		integer	:: i,j
		val=A(:,[(j,j=size(A,2),1,-1)])
	end function

	function kron(A,B) result(val)
		implicit none
		real	:: A(:,:),B(:,:)
		real	:: val(size(A,1)*size(B,1),size(A,2)*size(B,2))
		integer	:: i,j
		do i=1,size(A,2)
			do j=1,size(A,1)
				val((j-1)*size(B,1)+1:j*size(B,1),(i-1)*size(B,2)+1:i*size(B,2))=A(j,i)*B
			enddo
		enddo
	end function
	
	function trace(A) result(val)
		implicit none
		real	:: A(:,:)
		real	:: val
		integer	:: i
		val=0
		do i=1,size(A,1)
			val=val+A(i,i)
		enddo
	end function

	function vech(A) result(val)
		implicit none
		real	:: A(:,:)
		real	:: val(size(A,1)*(size(A,1)+1)/2)
		integer	:: i,j,c
		c=1
		do i=1,size(A,1)
			do j=i,size(A,1)
				val(c)=A(j,i)
				c=c+1
			enddo
		enddo
	end function
	
	function loadcsv(filename, headerflag) result(mat)
        character(len=*):: filename
		logical, optional	:: headerflag
        real, allocatable  :: mat(:,:)
        integer         :: ioerr,m,n,i,j,myunit,l1,l2
        character(:), allocatable  :: line
		real	:: nanx=0.0/0.0

!		nanx=0.0/0.0
 
        open (file=filename,status='old',iostat=ioerr,newunit=myunit)
        if(ioerr/=0) then
            print *,"loadcsv: file not found: ",filename
            stop
        endif
        m=0
        do
            read (myunit,*,iostat=ioerr)
            if (ioerr/=0) exit     
            m=m+1
        enddo
        rewind(myunit)
		call getline
        n=countitems(line)
        rewind(myunit)
		if(present(headerflag)) then
			if(headerflag) then
				m=m-1
				call getline
			endif
		endif
        allocate(mat(m,n))
        do i=1,m
			call getline
			l1=1
			do j=1,n
				l2=merge(index(line(l1:),","),len(line)-l1+2,j<n)
				if(l2==0) then 
					print *,"loadcsv: not enough fields in row",i
					stop
				endif
				if(l2==1) then 
					mat(i,j)=nanx
				else
					read(line(l1:l1+l2-2),*,iostat=ioerr) mat(i,j)
					if(ioerr>0) mat(i,j)=nanx
				endif
				l1=l1+l2
			enddo
		enddo
        close(myunit)
        
        contains

			function countitems(line) result(val)
				character(len=*)    :: line
				integer				:: val
				integer				:: i,n
				val=1
				n=len(line)
				do i=1,n
					if(line(i:i)==",") val=val+1
				enddo
			end function
			
			subroutine getline
		        character(len=8192) :: buffer
				integer				:: lsize
		        line=''
		        do
		            read(myunit,"(a)",advance="no",iostat=ioerr,size=lsize) buffer
				    if(ioerr>0) then
						print *,"loadmat can't read first line: ",filename
						stop
					endif
					line=line//buffer(:lsize)
					if(ioerr<0) exit
				enddo
			end subroutine
	end function


    function loadstrings(filename) result(vec)
        character(len=*):: filename
        character(len=:), allocatable  :: vec(:)
        integer         :: ioerr,m,n,i,lsize,myunit
        character(:), allocatable  :: line
        character(len=256)    :: buffer
 
        open (file=filename,status='old',iostat=ioerr,newunit=myunit)
        if(ioerr/=0) then
            print *,"loadstrings file not found: ",filename
            stop
        endif
        m=0; n=0
        do
            read (myunit,"(a)",iostat=ioerr) buffer
            if (ioerr/=0) exit     
            m=m+1
			n=max(n,len_trim(buffer))
        enddo
        rewind(myunit)
		allocate(character(n)::vec(m))
		do i=1,m
            read(myunit,"(a)",iostat=ioerr) vec(i)
            if (ioerr/=0) then
                print *,"loadmat trouble reading from file: ",filename
				print *,i
                stop
            endif
        enddo
        close(myunit)
    end function
	
    function loadvec(filename) result(vec)
        character(len=*):: filename
        real, allocatable  :: vec(:)
        vec=[loadmat(filename)]
    end function
	
    function loadmat(filename) result(mat)
        character(len=*):: filename
        real, allocatable  :: mat(:,:)
        integer         :: ioerr,m,n,i,lsize,myunit
        character(len=256)    :: buffer
        character(:), allocatable  :: line
 
        open (file=filename,status='old',iostat=ioerr,newunit=myunit)
        if(ioerr/=0) then
            print *,"loadmat file not found: ",filename
            stop
        endif
        m=0
        do
            read (myunit,*,iostat=ioerr)
            if (ioerr/=0) exit     
            m=m+1
        enddo
        rewind(myunit)
        line=''
        do
            read(myunit,"(a)",advance="no",iostat=ioerr,size=lsize) buffer
!             read(myunit,"(a)",advance="yes",iostat=ioerr,size=lsize) buffer
            if(ioerr>0) then
                print *,"loadmat can't read first line: ",filename
                stop
            endif
            line=line//buffer(:lsize)
            if(ioerr<0) exit
        enddo
        n=ntokens(line)
        rewind(myunit)
        allocate(mat(m,n))
        do i=1,m
            read(myunit,*,iostat=ioerr) mat(i,:)
            if (ioerr/=0) then
                print *,"loadmat trouble reading from file: ",filename
				print *,i
                stop
            endif
        enddo
        close(myunit)
        
        contains

        function ntokens(line) result(val)
            integer :: val
            character(len=*)    :: line
            integer             :: i, n, toks
            CHARACTER TAB, space
            TAB=CHAR(9)
			space=char(32)
    
            i = 1
            n = len_trim(line)
            toks = 0
            val = 0
            do while(i <= n)
               do while(line(i:i)==space.or.line(i:i) == ' '.or.line(i:i) == tab.or.line(i:i) == ',') 
                 i = i + 1
                 if (n < i) return
               enddo
               toks = toks + 1
               val = toks
               do
                 i = i + 1
                 if (n < i) return
                 if (line(i:i) == ' '.or.line(i:i) == tab.or.line(i:i) == ',') exit
               enddo
            enddo
        end function ntokens
    end function

    subroutine savemat(filename,mat)
        character(len=*):: filename
        real            :: mat(:,:)
        integer         :: ioerr,i,j,myunit
 
        open (file=filename,iostat=ioerr,newunit=myunit)
        if(ioerr/=0) then
            print *,"savemat could not open file for writing: ",filename
            stop
        endif
        
        do i=1,size(mat,1)
            do j=1,size(mat,2)
                write(myunit,'(ES27.18E3)',advance='no') mat(i,j)
            enddo
            write(myunit,*)
        enddo
        close(myunit)
    end subroutine
        
    subroutine savevec(filename,vec)
        character(len=*):: filename
        real            :: vec(:),mat(size(vec),1)
        mat(:,1)=vec
        call savemat(filename,mat)
	end subroutine
    
    subroutine savematcsv(filename,mat)
        character(len=*):: filename
        real            :: mat(:,:)
        integer         :: ioerr,i,j,myunit
 
        open (file=filename,iostat=ioerr,newunit=myunit)
        if(ioerr/=0) then
            print *,"savemat could not open file for writing: ",filename
            stop
        endif
        
        do i=1,size(mat,1)
            do j=1,size(mat,2)
                write(myunit,'(ES27.18E3)',advance='no') mat(i,j)
				if(j<size(mat,2)) write(myunit,'(a)',advance='no') ","
            enddo
            write(myunit,*)
        enddo
        close(myunit)
	end subroutine

	subroutine mydispnames(vec)
		character(len=*)	:: vec(:)
		character(len=size(vec)*10) :: line
		integer	:: i
		do i=1,size(vec)
			write(line((i-1)*10+1:i*10),'(a)') vec(i)
		enddo
		write (*,'(a)') line
	end subroutine
	
	subroutine mydispmat(mat)
		real	:: mat(:,:)
		integer	:: i
		
!$omp critical		
		do i=1,size(mat,dim=1)
			call mydispline(mat(i,:))
		enddo
		write (*,*)
!$omp end critical
	end subroutine

	subroutine mydispscalar(x)
		real	:: x
		write (*,'(a)') numstr(x)
		write (*,*) 
	end subroutine

	subroutine mydispstring(str)
		character(len=*)	:: str
		write (*,'(a)') str
	end subroutine

	
	subroutine mydispvec(vec)
		real	:: vec(:)
!$omp critical		
		call mydispline(vec)
		write (*,*) 
!$omp end critical
	end subroutine
	
	subroutine mydispline(vec)
		real	:: vec(:)
		integer	:: n,i
		character(len=size(vec)*10) :: line
		n=size(vec)
		do i=1,n
			line((i-1)*10+1:i*10)=numstr(vec(i))
		enddo
		write (*,'(a)') line
	end subroutine
	
	function numstr(x) result(val)
		real	:: x,ax
		character(len=10)	:: val
			
		ax=abs(x)
		if(ax==0) then
			write(val,'(I9)') 0
			return
		endif
		if(ax>0.001) then
			if(ax<9.999) then
				write(val,'(F9.6)') x
				return
			endif
			if(ax<99.999) then
				write(val,'(F9.5)') x
				return
			endif
			if(ax<999.99) then
				write(val,'(F9.4)') x
				return
			endif
			if(ax<9999.9) then
				write(val,'(F9.3)') x
				return
			endif
		endif
		if(ax<1) ax=1/ax
		if(ax<9.99E99) then
			write(val,'(ES9.2e2)') x
			return
		endif
		write(val,'(ES9.1e3)') x
	end function

	subroutine mydispivec(vec) 
		integer	:: vec(:)
		integer	:: n,i
		character(len=size(vec)) :: line
		
		n=size(vec)
		do i=1,n
			write(line(i:i),'(I1)') vec(i)
		enddo
		write (*,'(a)') line
	end subroutine
	
	function convtos(i,k) result(val)
		character(:),allocatable :: val
		integer	:: i
		integer, optional	:: k
		character(range(i)+2) :: tmp
		character(10)		:: fmt
		
		if(present(k)) then
			write (fmt,'(a,I0,a)') "(i0.",k,")"
			write (tmp,fmt) i
		else
		  write(tmp,'(i0)') i
		endif
		val = trim(tmp)
	end function
		
end module	
		
		
