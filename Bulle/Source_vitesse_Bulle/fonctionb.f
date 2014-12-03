c****************************************************************
c                fonction des termes compris dans la 
c                     Force a l'ordre 1 pour la
c                           sphere solide
c     
c
c     Appelée par la fonction gaussint1dim pour le calcul
c     de l'integrale de la force
c
c     Programmeur: M.Guemas; Novembre 2012.
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

      function fonctionb(eta,Nt,Nx,c,zetap,hxi,dhxi,tabxi,
     &             Bn,Dn,alpha,alpha0,b1,e1,d1,b2,e2,d2)
c     ==================================================================    
c     declaration des variables en entree
c     -----------------------------------
      integer, intent(in)  ::  Nt,Nx      
      double precision, intent(in)  :: eta, c, zetap 
      double precision :: alpha0
      double precision, dimension(Nt) :: alpha
      double precision, dimension(Nt+1) :: Bn, Dn
      double precision, dimension(Nx) :: hxi, dhxi
      double precision, dimension(0:Nx) :: tabxi   
      double precision, dimension(Nx) :: b1, e1, d1  
      double precision, dimension(Nx) :: b2, e2, d2     

c     declaration des variables locales
c     ---------------------------------
      integer :: n, k, i      
      double precision ::  dPn0, uEta, duZeta
      double precision ::  SigmaEta, SigmaZeta, dSigma
      double precision ::  dSigmazz
      double precision, dimension(Nt+3) :: Pn, dPn
      double precision, dimension(Nt+1) :: Un, dUn    
      double precision, dimension(Nt+1) :: d2Un, d3Un      
      double precision :: zeta1, dzeta1
      double precision :: res, dres, xk, xi,ds,dx, xn
      double precision :: hInt, dhInt   
   
c     Declaration des variables en sorties
c     ------------------------------------
      double precision :: fonctionb  

c     Declaration des interfaces
c     ---------------------------
      interface

         function u0etaB(Nt,c,zetap,xi,Pn,Un,dUn)
           integer, intent(in)  ::  Nt
           double precision, intent(in) :: xi, c, zetap           
           double precision, dimension(Nt+3), intent(in) :: Pn
           double precision, dimension(Nt+1), intent(in) :: Un
           double precision, dimension(Nt+1), intent(in) :: dUn  
           double precision :: u0etaB
         end function u0etaB

         function du0zetaB(Nt,c,zetap,xi,Pn,Un,dUn)
           integer, intent(in)  ::  Nt
           double precision, intent(in) :: xi, c, zetap
           double precision, dimension(Nt+3), intent(in) :: Pn  
           double precision, dimension(Nt+1), intent(in) :: Un
           double precision, dimension(Nt+1), intent(in) :: dUn  
           double precision :: du0zetaB
         end function du0zetaB

         function sigma0zzB(Nt,c,zetap,xi,alpha,alpha0,Pn,
     &                      Un,dUn)     
           integer, intent(in)  ::  Nt
           double precision, intent(in) :: c, zetap 
           double precision, intent(in) :: xi, alpha0      
           double precision, dimension(Nt), intent(in) :: alpha           
           double precision, dimension(Nt+3), intent(in) :: Pn        
           double precision, dimension(Nt+1), intent(in) :: Un
           double precision, dimension(Nt+1), intent(in) :: dUn  
           double precision :: sigma0zzB
         end function sigma0zzB

	 function sigma0xxB(Nt,c,zetap,xi,alpha,alpha0,Pn,
     &                      Un,dUn) 
           integer, intent(in)  ::  Nt
           double precision, intent(in) :: c, zetap 
           double precision, intent(in) :: xi, alpha0
           double precision, dimension(Nt), intent(in) :: alpha           
           double precision, dimension(Nt+3), intent(in) :: Pn        
           double precision, dimension(Nt+1), intent(in) :: Un  
           double precision, dimension(Nt+1), intent(in) :: dUn 
           double precision :: sigma0xxB
	 end function sigma0xxB

         function dsigma0zxB(Nt,c,zetap,xi,Pn,Un,dUn,d2Un,d3Un)
           integer, intent(in)  ::  Nt
           double precision, intent(in) :: xi, c, zetap                
           double precision, dimension(Nt+3), intent(in) :: Pn
           double precision, dimension(Nt+1), intent(in) :: Un
           double precision, dimension(Nt+1), intent(in) :: dUn
           double precision, dimension(Nt+1), intent(in) :: d2Un 
           double precision, dimension(Nt+1), intent(in) :: d3Un  
           double precision :: dsigma0zxB
	 end function dsigma0zxB

         function dsigma0zzB(Nt,c,zetap,xi,Pn,Un,dUn,d2Un)
           integer, intent(in)  ::  Nt
           double precision, intent(in) :: xi, c, zetap                
           double precision, dimension(Nt+3), intent(in) :: Pn
           double precision, dimension(Nt+1), intent(in) :: Un
           double precision, dimension(Nt+1), intent(in) :: dUn
           double precision, dimension(Nt+1), intent(in) :: d2Un             
           double precision :: dsigma0zzB
	 end function dsigma0zzB

	 
         function ispline(u, x, y, b, e, d, n)
           integer n
           double precision  u,x(n),y(n),b(n),e(n),d(n)
           double precision ispline
         end function ispline

      end interface
c     ==================================================================      

      
c     On appel le sous programme des Pn,Bn,Dn  au debut
c     du code car ils seront utilisé plusieurs fois

c     Soit xi en fonction de l'angle eta entre 0 et pi
c     -----------------------------------------------
       xi=dcos(eta) 
c      write(*,*) 'xi=',xi
       
       if (eta.eq.0) then
            open(unit=8,file='bloc31.dat',status='replace')
         else
            open(unit=8,file='bloc31.dat',position='append')
         endif 
         
c     On initialise Pn, dPn, dUn et d3Un à O
c     --------------------------------------
      do 1 k=1,Nt+3
       Pn(k)=0.d0
       dPn(k)=0.d0
 1    continue

      do 2 k=1,Nt+1  
       Un(k)=0.d0
       dUn(k)=0.d0
       d2Un(k)=0.d0
       d3Un(k)=0.d0       
 2    continue

c     Construction du tableau Pn,dPn,dUn et d3Un 
c     ------------------------------------------
c     On sait que Pn0=1.D0      
      call DLegendre(0,xi,dres)      
      dPn0=dres
      do 9 n=1,Nt+3 
        xn=dble(n)       
        call Legendre(n,xi,res)
        call DLegendre(n,xi,dres)
        Pn(n)=res
        dPn(n)=dres
        write(8,*) n, Pn(n)
 9      continue     
    
      do 10 n=1,Nt+1 
        xn=dble(n)      
        Un(n)=Bn(n)*dsinh((xn-5.d-1)*zetap)
     &     +Dn(n)*dsinh((xn+(3.d0/2.d0))*zetap)
        dUn(n)=(xn-5.d-1)*Bn(n)*dcosh((xn-5.d-1)*zetap)
     &         +(xn+(3.d0/2.d0))*Dn(n)
     &         *dcosh((xn+(3.d0/2.d0))*zetap) 
        d2Un(n)=((xn-(1.d0/2.d0))**2.d0)*Bn(n)
     &          *dsinh((xn-(1.d0/2.d0))*zetap)
     &          +((xn+(3.d0/2.d0))**2.d0)*Dn(n)
     &          *dsinh((xn+(3.d0/2.d0))*zetap) 
        d3Un(n)=((xn-5.d-1)**3.d0)*Bn(n)
     &          *dcosh((xn-5.d-1)*zetap)
     &          +((xn+(3.d0/2.d0))**3.d0)*Dn(n)
     &          *dcosh((xn+(3.d0/2.d0))*zetap) 
c        write(8,*) n, Bn(n),Dn(n) 
c        write(8,*)  n, Un(n), dUn(n),d2Un(n),d3Un(n)
c        write(8,*)'d2Un(',n,')',d2Un(n)
 10   continue  


c$$$      do n=1,Nt
c$$$         write(16,*) n,Bn(n),Dn(n),Un(n),dUn(n),alpha(n),
c$$$     &               alpha0,Pn(n)
c$$$      enddo

c2345678912345678912345678912345678912345678912345678912345678912
c	 On interpole h(xi) et h'(x)
c	 ---------------------------
c$$$        open(unit=7,file='coef.dat')
c$$$        write(7,*)'Coefficients 1 :'
c$$$        write(7,*)        
c$$$         do 111 i=1,Nx  
c$$$c         write(*,*) i, b1,e1,d1
c$$$         write(7,*) i, b1(i), e1(i), d1(i)
c$$$111       continue
c$$$        write(7,*)
c$$$        write(7,*)'Coefficients 2 :'
c$$$        write(7,*)
c$$$         do 110 i=1,Nx           
c$$$         write(7,*) i, b2(i), e2(i), d2(i)
c$$$ 110     continue   


      ! ATTENTION on doit refaire le tableau tabxi car il a été ecrasé par la variable
      ! x(i) dans la fonction gaussintdim.f

c$$$      tabxi(0)=-1.d0
c$$$      do 12 n=0,Nx
c$$$        tabxi(i)=-1.D0+dx*dble(i)
c$$$ 12   continue   
c$$$      tabxi(Nx)=1.d0
   
	 hInt=ispline(xi,tabxi,hxi,b1,e1,d1,Nx)           
	 dhInt=ispline(xi,tabxi,dhxi,b2,e2,d2,Nx)
c         write(*,*) 'hInt=', hInt
c         write(*,*) 'dhInt=', dhInt
    
c         write(8,*) 'uEta_barre=', uEta_barre
c         write(8,*) eta,dSigma_barre
c         write(8,*)'SigmaEta_barre=', SigmaEta_barre
c         write(8,*) 'SigmaZeta_barre=', SigmaZeta_barre
c         write(8,*) 'duZeta_barre=', duZeta_barre
c         write(8,*) 'dSigma_barre=', dSigma_barre

c     Variables contenant le resultats des fonctions vitesse et contrainte
c     --------------------------------------------------------------------
c     u0etaSL
c         write(*,*)'c avant ueta',c
      uEta=u0etaB(Nt,c,zetap,xi,Pn,Un,dUn)
c     du0zetaSL
      duZeta=du0zetaB(Nt,c,zetap,xi,Pn,Un,dUn)
c     sigma0xxSL
      SigmaEta=sigma0xxB(Nt,c,zetap,xi,alpha,alpha0,
     &                   Pn,Un,dUn)
c     sigma0zzSL
      SigmaZeta=sigma0zzB(Nt,c,zetap,xi,alpha,
     &                    alpha0,Pn,Un,dUn)
c     dsigma0zx
      dSigma=dsigma0zxB(Nt,c,zetap,xi,Pn,Un,dUn,d2Un,d3Un)
      
      dSigmazz=dsigma0zzB(Nt,c,zetap,xi,Pn,Un,dUn,d2Un) 

c      write(8,*) xi, uEta
c          write(8,*) xi, SigmaZeta
c     &    (2.d0*SigmaZeta_barre-SigmaEta_barre)
c     &    *(-(1.d0-xi)*dhInt+(3.d0/2.d0)*hInt) 
c          write(8,*) eta,
c     &         (-(1.d0-xi)*dhInt+(3.d0/2.d0)*hInt)
c         write(8,*) eta,
c     &       -hInt*duZeta_barre*SigmaZeta_barre 
c           write(8,*) eta,((1.d0-xi)**2)*hInt
c     &       *uEta_barre*dSigma_barre
c$$$     &               

c
c     Soit l'expression de la fonctionb
c       fonctionb = {u0etaSL*[(1-x)^(3/2)*dhInt-(3/2)*(1-x)^(1/2)*hInt]*[2*simga0zzSL-sigma0xxSL]
c                     +(1-x)^(3/2)*hInt*[u0etaSL*dsigma0z - sigma0zz*du0zeta]} dS

c     Expression de zeta1 et dzeta1
c     -----------------------------
      zeta1=((1.d0-xi)**(3.d0/2.d0))*hInt
      dzeta1=((1.d0-xi)**(3.d0/2.d0))*dhInt
     &       -(3.d0/2.d0)*((1.d0-xi)**(1.d0/2.d0))*hInt
      ds=dsin(eta)/((1.d0-xi)**2.d0)

c      Expression de la fonction fb
c      ----------------------------
c$$$      fonctionb=(uEta*dzeta1*(2.d0*SigmaZeta
c$$$     &           -SigmaEta)+zeta1*(uEta*
c$$$     &           dSigma-SigmaZeta*duZeta))*ds           
      
      !test
      fonctionb=1.d-3
      close(8)  
           
      end function fonctionb
c     -------- Fin du programme----------------------------------
c****************************************************************

