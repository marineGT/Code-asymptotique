c****************************************************************
c                fonction des termes compris dans la 
c                     Force a l'ordre 1 pour la
c                           sphere solide
c     
c
c     AppelÃ©e par la fonction gaussint1dim pour le calcul
c     de l'integrale de la force
c
c     Programmeur: M.Guemas; Novembre 2012.
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

      function verifsigma0zzB(teta,Nt,c,zetap,alpha,alpha0,
     &                        Bn,Dn)
c     ==================================================================    
c     declaration des variables en entree
c     -----------------------------------
      integer, intent(in)  ::  Nt      
      double precision, intent(in)  :: teta
      double precision, intent(in) ::  c, zetap
      double precision, dimension(Nt) :: alpha
      double precision, intent(in) :: alpha0 
      double precision, dimension(Nt+1) :: Bn,Dn

c     declaration des variables locales
c     ---------------------------------
      integer :: n, k
      double precision, dimension(Nt+3) :: Pn
      double precision, dimension(Nt+1) :: Un, dUn
      double precision :: res, xn, ch, xi
   
c     Declaration des variables en sorties
c     ------------------------------------
      double precision :: verifsigma0zzB

c     Declaration des interfaces
c     ---------------------------
      interface

         function sigma0zzB(Nt,c,zetap,xi,alpha,alpha0,Pn,
     &                      Un,dUn)     
           integer, intent(in)  ::  Nt
           double precision, intent(in) :: c, zetap 
           double precision, intent(in) :: xi
           double precision, dimension(Nt), intent(in) :: alpha 
           double precision, intent(in) :: alpha0 
           double precision, dimension(Nt+3), intent(in) :: Pn        
           double precision, dimension(Nt+1), intent(in) :: Un
           double precision, dimension(Nt+1), intent(in) :: dUn  
           double precision :: sigma0zzB
         end function sigma0zzB

        subroutine Legendre(n,x,resultat)
           integer :: n
           double precision :: x,resultat           
       end subroutine Legendre

      end interface
c     ==================================================================      
     
c     On appel le sous programme des Pn,Bn,Dn  au debut
c     du code car ils seront utilisé plusieurs fois

c     Soit xi en fonction de l'angle eta entre 0 et pi
c     -----------------------------------------------
     
      ch=dcosh(zetap)
      xi=(1.d0-dcos(teta)*ch)/(ch-dcos(teta)) 
c      write(*,*) 'teta=', teta, 'xi=', xi 

c     On initialise Pn, Un et dUn 
c     --------------------------------------
      do 1 k=1,Nt+3
       Pn(k)=0.d0
 1    continue

      do 2 k=1,Nt+1  
       Un(k)=0.d0
       dUn(k)=0.d0
 2    continue

c     Construction du tableau Pn,dPn,dUn et d3Un 
c     ------------------------------------------
      do n=1,Nt+3 
        xn=dble(n)       
        call Legendre(n,xi,res)
        Pn(n)=res
c        write(17,*) xi, n, Pn(n)
      enddo
    
      do 10 n=1,Nt+1 
        xn=dble(n)      
        Un(n)=Bn(n)*dsinh(-(xn-5.d-1)*zetap)
     &     +Dn(n)*dsinh(-(xn+(3.d0/2.d0))*zetap)
        dUn(n)=(xn-5.d-1)*Bn(n)*dcosh((xn-5.d-1)*zetap)
     &         +(xn+(3.d0/2.d0))*Dn(n)
     &         *dcosh((xn+(3.d0/2.d0))*zetap) 
c        write(16,*) xi, n, Un(n), dUn(n)
 10   continue  


c     Calcul de sigma0zz pour une valeur de teta
c     ------------------------------------------

      verifsigma0zzB=sigma0zzB(Nt,c,zetap,xi,alpha,
     &               alpha0,Pn,Un,dUn)
      verifsigma0zzB=verifsigma0zzB*dsin(teta)*dcos(teta)

c      write(16,*) teta, xi, verifsigma0zzB

      end function verifsigma0zzB
c     -------- Fin du programme----------------------------------
c****************************************************************

