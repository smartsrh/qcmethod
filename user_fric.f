!-------------------------------------------------------------------------------
!
!  Quasicontinuum (QC) Method: Mixed Continuum and Atomistic Simulation Package
!  QC Package distribution version 1.4 (November 2011)
!
!  Copyright (C) 2003 R. Miller, M. Ortiz, R. Phillips, D. Rodney, E. B. Tadmor
!  Copyright (C) 2004, 2005, 2006, 2007, 2011, 2011 R. Miller, E. B. Tadmor
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
!
!  For more information please visit the QC website at:
!
!    www.qcmethod.com
!
!  or contact Ron Miller or Ellad Tadmor with contact information below:
!
!  NOTE HOWEVER THAT NEITHER RON MILLER NOR ELLAD TADMOR WILL BE ABLE TO
!  PROVIDE SUPPORT OR DEBUGGING SERVICES FOR THE CODE. THE CODE IS SUPPLIED 
!  AS IS FOR YOUR USE.
!
!    Ron Miller
!    Department of Mechanical and Aerospace Engineering
!    Carleton University
!    1125 Colonel By Drive
!    Ottawa, ON, K1S 5B6
!    CANADA
!   
!    URL: http://www.mae.carleton.ca/Ron_Miller
!
!  -or-
!
!    Ellad B. Tadmor
!    Department of Aerospace Engineering and Mechanics
!    University of Minnesota
!    107 Akerman Hall, 110 Union St SE
!    Minneapolis, MN 55455
!    USA
!    
!    URL: http://www.aem.umn.edu/people/faculty/bio/Tadmor.shtml
!   
!-------------------------------------------------------------------------------

!****************************************************************************
!**
!**  MODULE mod_uservars   contains variables used by user routines
!**
!**                        Friction example
!**
!**  Variable Definitions:
!**  ---------------------
!**
!**  Integer Variables:
!**
!**    user_iparam          -  array of integers available for storage
!**                            of user-specified parameters. See also
!**                            user_rparam.
!**
!**  Real Variables:
!**
!**    user_rparam          -  array of reals available for storage
!**                            of user-specified parameters. See also
!**                            user_iparam.
!**
!****************************************************************************

      module mod_uservars
      use mod_global, only : dp

      save
      private
      public user_iparam, user_rparam

      integer user_iparam(10)
      real(kind=dp) user_rparam(10)

      end module mod_uservars


!**-----------------------------------------------------------------------------
!**
!**   user_mesh : Mesh generator routine for friction example
!**
!**   Calling format:
!**
!**   mesh,,nx,ny
!**
!**   where  nx = number of x-divisions
!**          ny = number of y-divisions
!--
      subroutine user_mesh(id,x,ix,f,b,itx,key,input)
      use mod_uservars
      use mod_global,   only : dp,iregion,maxel,maxnp,ndf,ndm,nen,
     &                         nregion,numel,numnp,nxdm,SymmetricMesh
      use mod_output,   only : erroroutput,nerror,output,string,error,
     &                         logoutput
      use mod_boundary, only : ncb,nce,NCEMAX,elist
      use mod_repatom,  only : protected
      use mod_qclib,    only : next,freein
      use mod_mesh,     only : delaunay,PerturbMesh
      use mod_grain,    only : GetGrainNumVrts,GetGrainVertex,
     &                         NearestBSite
      implicit none

!-- Transferred variables
      integer,           intent(inout)  :: id(ndf,maxnp),
     &                                     ix(nen,maxel),
     &                                     itx(3,maxel)
      real(kind=dp),     intent(inout)  :: x(nxdm,maxnp),
     &                                     f(ndf,maxnp),
     &                                     b(ndf,maxnp)
      character(len=4),  intent(in)     :: key
      character(len=80), intent(in)     :: input

!-- Local variables
      real(kind=dp), allocatable :: vtmp(:)
      logical xflag,yflag,PointInPoly
      integer nx,ny,kount,igr,i,j,k,idum,kk,ivert,lower,upper,j2,
     &     numvrts,endofout
      real(kind=dp) xlo,xhi,ylo,yhi,dx,dy,xx,yy,xtmp(3),xtmp1(2),
     &     xmin,xmax,ymin,ymax,xc,yc,small,dum,vec(2),vlen,vecnorm(2)
     $     ,frac,xtri(2,3),aspect,vertex(ndm),vertex_next(ndm)
      real(kind=dp), parameter :: tol=1.e-3_dp

c     Parse input
      lower = 0
      upper = next(lower,input)
      call freein(input,lower,upper,nx,dum,1)
      lower = upper
      upper = next(lower,input)
      call freein(input,lower,upper,ny,dum,1)


      string=''; call output
      write(string,'(a)')
     &     'Generating uniform mesh for friction test problem'
      call output
      string=''; call output
      write(string,'(a,i4)')   '# x-elements           = ',nx
      call output
      write(string,'(a,i4)')   '# y-elements           = ',ny
      call output

      nregion=1
      numnp=0
      do i=1,1
         idum=numnp+1

         call GetGrainNumVrts(i,numvrts)

C        Dealing with the outer first iterate each vertex, 
C        get outer boundary nodes with linear interpolation 
C        in counter-clockwise

         do j=1,numvrts/2
            j2=j+1
            if(j2.gt.(numvrts/2)) j2=1
            call GetGrainVertex(i,j,vertex)
            call GetGrainVertex(i,j2,vertex_next)
            vec=vertex_next(1:2)-vertex(1:2)
            vlen=sqrt(dot_product(vec,vec))
            vecnorm=vec/vlen
            frac=vlen/10.0_dp

            string=''; call output
            string = 'display vec:'; call output
            string = ' from    to    x     y'
            call output
            write(string,'(i5,2x,i5,2x,f12.5,2x,f12.5)')j,j2,vec(1:2)
            call output

            do k=0,9
               xtmp1=vertex(1:2)+k*frac*vecnorm
               call NearestBSite(xtmp1,.false.,xtmp,igr,idum,x,b,ix,itx)
               do kk=1,numnp
                  vec=x(1:2,kk)-xtmp(1:2)
                  if(dot_product(vec,vec).lt.0.01_dp) go to 1020
               enddo

               string=''; call output
               string = 'display nodes on the last edge:'; call output
               write(string,'(i5,2x,f12.5,2x,f12.5)')numnp+1,xtmp(1:2)
               call output

               numnp = numnp + 1
               iregion(numnp)=igr
               x(1:3,numnp) = xtmp(1:3)
               nce=nce+1
               if(nce.gt.NCEMAX) stop 'stop: NCEMAX'
               elist(1,nce)=numnp
               elist(2,nce)=numnp+1
 1020          continue
            enddo          
         enddo
         elist(2,nce)=idum
         endofout=nce
         
C        then with inner, 
C        iterate each vertex, 
C        get inner boundary nodes with linear interpolation in clockwise, 

         do j=(numvrts)/2+1, numvrts
            j2=j+1
            if(j2.gt.(numvrts)) j2 = numvrts / 2 + 1
            call GetGrainVertex(i,j,vertex)
            call GetGrainVertex(i,j2,vertex_next)
            vec=vertex_next(1:2)-vertex(1:2)
            vlen=sqrt(dot_product(vec,vec))
            vecnorm=vec/vlen
            frac=vlen/10.0_dp
            do k=0,9
               xtmp1=vertex(1:2)+k*frac*vecnorm
               call NearestBSite(xtmp1,.false.,xtmp,igr,idum,x,b,ix,itx)
               do kk=1,numnp
                  vec=x(1:2,kk)-xtmp(1:2)
                  if(dot_product(vec,vec).lt.0.01_dp) go to 1030
               enddo
               numnp = numnp + 1
               iregion(numnp)=igr
               x(1:3,numnp) = xtmp(1:3)
               nce=nce+1
               elist(1,nce)=numnp
               elist(2,nce)=numnp+1
 1030          continue
            enddo            
         enddo

C        and make the last node in `elist` 
C        point to the start point of the iteration

         elist(2,nce) = endofout + 1

      enddo

      string=''; call output
      string = 'elist List:'; call output
      string = ' elist    start  end'
      call output
      do k=1,nce
        write(string,'(i5,2x,i5,2x,i5,2x,f12.5,2x,f12.5)')k,elist(1:2,k)
        call output
      enddo
      string = 'Node List:'; call output
      string = ' node        x             y     '
      call output
      do k=1,numnp
        write(string,'(i5,2x,f12.5,2x,f12.5)')k,x(1:2,k)
        call output
      enddo

      ncb=nce
      kount=numnp
c
c generate internal nodes
c
      xlo=1.e30_dp
      xhi=-1.e30_dp
      ylo=1.e30_dp
      yhi=-1.e30_dp
      do igr=1,1
         call GetGrainNumVrts(igr,numvrts)
         do ivert=1,numvrts
            call GetGrainVertex(igr,ivert,vertex)
            xx = vertex(1)
            yy = vertex(2)
            xlo = min(xx,xlo)
            xhi = max(xx,xhi)
            ylo = min(yy,ylo)
            yhi = max(yy,yhi)
         enddo
      enddo

      dx = (xhi - xlo)/real(nx,dp)
      dy = (yhi - ylo)/real(ny,dp)

c     Generate Nodes
      do i = 1,nx-1
      do j = 1,ny-1
         xtmp1(1) = xlo + real(i,dp)*dx
         xtmp1(2) = ylo + real(j,dp)*dy
         if (abs(xtmp1(2))<0.1) goto 1010
         call NearestBSite(xtmp1,.false.,xtmp,igr,idum,x,b,ix,itx)
C          do k=1,kount
C             vec=x(1:2,k)-xtmp(1:2)
C             if (dot_product(vec,vec)<0.01) goto 1010
C          enddo         
         numnp = numnp + 1
         iregion(numnp)=igr
         x(1:3,numnp) = xtmp(1:3)
 1010    continue
      enddo
      enddo

c     Identify model boundaries
      xmin = x(1,1)
      xmax = x(1,1)
      ymin = x(2,1)
      ymax = x(2,1)
      do i = 1,numnp
         xx = x(1,i)
         yy = x(2,i)
         if (xx.lt.xmin) xmin=xx
         if (xx.gt.xmax) xmax=xx
         if (yy.lt.ymin) ymin=yy
         if (yy.gt.ymax) ymax=yy
      enddo
      xc = (xmin+xmax)/2.0_dp
      yc = (ymin+ymax)/2.0_dp

      ! Check for repeating nodes
      do i = 1, numnp
         do j = i+1, numnp
            if ( abs( x(1,i) - x(1,j)) < tol .and.
     &           abs( x(2,i) - x(2,j)) < tol ) then
               nerror = 2
               write(error(1),'(a,i5,a,i5)')'Node',j,' coincides with',i
               error(2) =
     &             'Node list sent to log file (subroutine user_mesh)'
               call erroroutput
               ! list nodes to log file
               string = 'Node List:'; call logoutput
               string = ' node        x             y     '
               call logoutput
               do k=1,numnp
                  write(string,'(i5,2x,f12.5,2x,f12.5)')k,x(1:2,k)
                  call logoutput
               enddo
               stop
            endif
         enddo
      enddo

c     Triangulate
      string=''; call output
      write(string,'(a,i5)')'**numnp = ',numnp
      call output

      SymmetricMesh = .false.   ! turn on symmetric meshing option


      call PerturbMesh(x,.true.)
      call delaunay(x,ix,b,f,id,itx)
      call PerturbMesh(x,.false.)
c
c delete elements with a poor aspect ratio
c
c      do i=numel,1,-1
c         do j=1,3
c            xtri(1:2,j)=x(1:2,ix(j,i))
c         enddo
c         call getaspect(xtri,aspect)
c         if(aspect.lt.1.e-1) then
c            ix(1:nen,i)=ix(1:nen,numel)
c            numel=numel-1
c         endif
c      enddo
c      
c      do i=1,numel
c         call CreateItx(i,ix,itx)
c      enddo

c      
c     Boundary conditions
      ymax=ymax-0.1_dp
      ymin=ymin+0.1_dp
      xmax=xmax-0.1_dp
      xmin=xmin+0.1_dp

      do i = 1,numnp
         xx=x(1,i)
         yy=x(2,i)
c        Protect initial nodes from deletion during coarsening
         protected(i) = .true.
c        Rigid boundary at x=xmax
          if (xx.gt.xmax) then
             id(1,i)=1
C               f(1,i)=1.0_dp
             id(2,i)=1
             id(3,i)=1
          endif
C         Rigid boundary at x=xmin
          if (xx.lt.xmin) then
             id(1,i)=1
             id(2,i)=1
             id(3,i)=1
          endif
      enddo

c     Report and exit
      write(string,'(a,i5)')'**numel = ',numel
      call output

      user_rparam(4) = xmax
      user_rparam(5) = xmin

      return
      end


!**---------------------------------------------------------------
!**   user_bcon : incremental uniform shear in x-direction
!**
!**   Calling format:
!**
!**        bcon
!**
!--
      subroutine user_bcon(id,x,ix,f,jdiag,str,eps,q,b,dr,db,shp,
     &   xsj,key,input,flag)
      use mod_uservars
      use mod_global,  only : dp,maxel,maxneq,maxnp,ndf,nen,nq,nstr,
     &                        numnp,nxdm
      use mod_output,  only : output,string
      use mod_pload,   only : GetLoadPropFact,GetLoadOldPropFact,
     &                        GetLoadTime,GetLoadTimeStep
      use mod_grain,   only : PointInGrain
      implicit  none

!-- Transferred variables
      real(kind=dp),     intent(inout) :: b(ndf,maxnp),str(nstr,maxel),
     &                                    eps(nstr,maxel),q(nq,maxel),
     &                                    x(nxdm,maxnp),f(ndf,maxnp),
     &                                    dr(ndf,maxnp),db(maxneq),
     &                                    shp(3,nen,maxel),xsj(maxel)
      integer,           intent(inout) :: jdiag(maxneq),ix(nen,maxel),
     &                                    id(ndf,maxnp)
      character(len=80), intent(in)    :: input
      character(len=4),  intent(in)    :: key
      logical,           intent(inout) :: flag

!-- Local variables
C       integer i
C       real(kind=dp) delta,propfact,propolfact

C       ! proportional load factors
C       propfact = GetLoadPropFact()
C       propolfact = GetLoadOldPropFact()
C       delta=propfact-propolfact

C       do i=1,numnp
C          if(PointInGrain(x(1:2,i),2)) then
C             b(1,i)=b(1,i)+delta
C          endif
C       enddo

      real(kind=dp), parameter :: tole=1.e-7_dp
      real(kind=dp) xmax,xmin
      real(kind=dp) xx,yy,propfact,delta,propolfact
      integer i

      ! Extract punch information stored in user param array 
      ! by meshing routine
      xmax = user_rparam(4)
      xmin = user_rparam(5)


      ! Constrain displacements of all nodes under indenter to
      ! appropriate values. List of nodes in contact with the indenter
      ! may change due to adaption - so need to look for them.
      propfact = GetLoadPropFact()
      propolfact = GetLoadOldPropFact()
      delta=propfact-propolfact
      do i=1,numnp
         if(PointInGrain(x(1:2,i),1)) then
           xx=x(1,i)
           yy=x(2,i)
           b(1,i)=b(1,i)+delta
           if (xx.gt.xmax) then
              id(1,i)=1
              id(2,i)=1
              id(3,i)=1
              f(1,i)=1.0_dp
              b(1,i)=propfact
           endif
           if (xx.le.xmin) then
              id(1,i)=1
              id(2,i)=1
              id(3,i)=1
              f(1,i)=-1.0_dp
              b(1,i)=-propfact
           endif
         endif
      enddo

      ! output
      write(string,'(a)')'* Boundary conditions applied:'
      call output
      write(string,'(2x,a,f5.2)')
     &             'nodes under indenter displaced down by ',propfact
      call output

      end


!**---------------------------------------------------------------
!**   user_pdel : compute the total external load for the p-delta 
!**               curve and for convergence checks.
!**
!--
      subroutine user_pdel(id,x,ix,f,jdiag,str,eps,q,b,dr,db,shp,
     &   xsj,key,input,force)
      use mod_uservars
      use mod_global, only : dp,maxel,maxneq,maxnp,ndf,nen,nq,nstr,
     &                       numnp,nxdm
      implicit  none

!-- Transferred Variables
      real(kind=dp),     intent(inout) :: b(ndf,maxnp),str(nstr,maxel),
     &                                    eps(nstr,maxel),q(nq,maxel),
     &                                    x(nxdm,maxnp),f(ndf,maxnp),
     &                                    dr(ndf,maxnp),db(maxneq),
     &                                    shp(3,nen,maxel),xsj(maxel)
      integer,           intent(inout) :: jdiag(maxneq),ix(nen,maxel),
     &                                    id(ndf,maxnp)
      character(len=80), intent(in)    :: input
      character(len=4),  intent(in)    :: key
      real(kind=dp),     intent(out)   :: force

!-- Local Variables
      integer i
      real(kind=dp) ff
 
      force = 0.0_dp
      do i = 1,numnp
         if (id(1,i).eq.1.and.x(2,i).gt.0.0_dp) force=force-dr(1,i)
      enddo

      return
      end

!**---------------------------------------------------------------
!**   user_potential : user-defined external potential 
!**
!**
!**
      subroutine user_potential(id,x,ix,f,jdiag,str,eps,q,b,dr,shp,xsj,
     &                          fls,flt,flw)
      use mod_uservars
      use mod_global, only : dp,maxneq,maxnp,maxel,ndf,nen,nq,nstr,nxdm,
     &                       user_energy
      use mod_pload,  only : GetLoadPropFact,GetLoadOldPropFact,
     &                       GetLoadTime,GetLoadTimeStep
      use mod_stiff,  only : tang,istiff
      implicit  none

!-- Transferred Variables
      real(kind=dp), intent(inout) :: b(ndf,maxnp),str(nstr,maxel),
     &                                eps(nstr,maxel),q(nq,maxel),
     &                                x(nxdm,maxnp),f(ndf,maxnp),
     &                                dr(ndf,maxnp),shp(3,nen,maxel),
     &                                xsj(maxel)
      integer,       intent(inout) :: jdiag(maxneq),ix(nen,maxel),
     &                                id(ndf,maxnp)
      logical,       intent(in)    :: fls,flt,flw

!-- Local Variables

      user_energy = 0.0_dp

      return
      end

!--------------------------------------------------------------------
!
!     user_plot : User-specified output plot
!
!     Calling format:
!
!     plot,user,filepref,[index],[umag],[scale],[append]
!
!
!
      subroutine user_plot(id,x,ix,itx,f,jdiag,str,eps,q,b,dr,db,shp,
     &                     xsj,index,umag,scale,logic,input,upper)
      use mod_uservars
      use mod_global, only : dp,maxnp,maxel,maxneq,ndf,nen,nq,nstr,nxdm
      implicit  none

!-- Transferred Variables
      real(kind=dp),     intent(in)    :: b(ndf,maxnp),str(nstr,maxel),
     &                                    eps(nstr,maxel),q(nq,maxel),
     &                                    x(nxdm,maxnp),f(ndf,maxnp),
     &                                    dr(ndf,maxnp),db(maxneq),
     &                                    shp(3,nen,maxel),xsj(maxel),
     &                                    umag,scale
      integer,           intent(in)    :: jdiag(maxneq),ix(nen,maxel),
     &                                    id(ndf,maxnp),itx(3,maxel),
     &                                    index,logic
      integer,           intent(inout) :: upper
      character(len=80), intent(in)    :: input

      return
      end
