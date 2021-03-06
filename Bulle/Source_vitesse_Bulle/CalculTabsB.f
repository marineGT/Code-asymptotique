c****************************************************************
c         Programme de calcul en coordonnes bipolaire des 
c       coefficients S_n dans le cas Bond=1 et Ca est petit
c         pour une bulle en translation perpendiculaire a une 
c           paroi plane avec resolution d'un systeme 
c           
c
c     Conventions:
c     La paroi plane est en z=0 et contient l'origine: c'est donc
c     le plan (O,x,y).
c     alti: "altitude" du centre de la particule par rapport a la
c           paroi solide
c     Nt: entier qui fixe la troncature dans la resolution 
c         de collocation en coordonnees multipolaires avec
c     Bond: nombre de Bond
c     lsa : altitude centre de la sphere/rayon a de la sphère
c     Nt: indice de troncature 
c     coefA: coefficiant A appartenant à la solution generale de 
c            l'equation differentielle voir fonction zeta 
c     zeta : fonction liee a la deformation de l'interface
c     eta : coordonnee bipsherique (bipolaire), est un angle compris 
c           entre 0 et pi      
c     xi : correspond a cos(eta) 
c     c : parametre liee a la distance l entre l'interface et la 
c         sphere, ie sqrt(l²-1)
c     Bn,Dn : coefficient lies a la fonction Un
c     zetap : indique que l'on se situe a la surface de la sphere 
c     kn: facteur utilise dans la definition des Bn et Dn
c     alpha: coefficient liee au la pression, utilise ici pour 
c            calculer les S_n
c     alpha0 : terme de alpha lorsque la pression est nulle a  
c              l'infini
c     bn : inconnus auxiliaires definis par les Bn et Dn servant
c          a calculer les alpha
c     dUn: derivee de la fonction Un(zeta) liee a la fonction de 
c          Stokes Psi
c
c     tabs(n=1,...,Nt): coefficients S_n de la solution particulière 
c                       de la deformee du problem
c
c     Programmeur: M.Guemas; Decembre 2011.
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

      subroutine CalculTabsB(lsa,Nt,tabs,tabs0,tabsB,tabsB0,Bn,
     &                       Dn,alpha,alpha0,R0,Rn)

c     ===========================================================
      implicit none 

c     declaration des variables en entrees
c     ------------------------------------
      integer, intent(in) :: Nt
      double precision, intent(in) :: lsa

c     Declaration des variables locales
c     ---------------------------------
      integer :: n,m
      double precision :: xn, xNt, kn, c, zetap
      double precision :: sum, xm, ch, sh, an 
      double precision,dimension(Nt+1) :: dUn, dUnB
      double precision,dimension(Nt) :: tabbn
      double precision,dimension(Nt) :: taban

c     declaration des variables en sorties
c     ------------------------------------
      double precision, intent(out) :: alpha0, tabs0
      double precision, intent(out) :: tabsB0, R0       
      double precision,dimension(Nt+1) :: Bn,Dn,tabs,tabsB
      double precision,dimension(Nt)  :: alpha,Rn
     
     
c     initialisation des tableaux Dn,Bn,alpha,alpha0,dU 
c     -------------------------------------------------
c$$$      do 1 k=1,3000
c$$$       tabs(k)=0.d0
c$$$ 1    continue
c      write(*,*)'On rentre dans CalculTabs'
c     calcul du reel c
c     ----------------
      c=dsqrt(lsa**2-1.d0)
c      write(*,*)'Nt en  entree  = ', Nt
c     calcul du reel zetap
c     --------------------
      zetap=-dlog(lsa+dsqrt(lsa**2-1.d0)) 
      ch=dcosh(zetap)
      sh=dsinh(zetap)
c     ___________________________________________________________
c     ------- Calcul des coefficients intermediaires ------------

c     construction du tableau Bn,Dn,kn et dUn pour n=1,Nt+1
c     -----------------------------------------------------
c$$$      do 2 n=1,Nt+1
c$$$        xn=dble(n)         
c$$$        an=(c**2)*xn*(xn+1.d0)/(sqrt(2.d0)*(2.d0*xn+1.d0))
c$$$        Dn(n)=-an*dexp(-(xn+(3.d0/2.d0))*zetap)
c$$$     +        /((2.d0*xn+3.d0)*dsinh((xn+(3.d0/2.d0))*zetap))     
c$$$        Bn(n)=an*dexp(-(xn-(1.d0/2.d0))*zetap)
c$$$     +        /((2.d0*xn-1.d0)*dsinh((xn-(1.d0/2.d0))*zetap))        
c$$$        Un(n)=-Bn(n)*dsinh((xn-(1.d0/2.d0))*zetap)-Dn(n)
c$$$     &         *dsinh((xn+(3.d0/2.d0))*zetap) 
c$$$        dUn(n)=(xn-5.d-1)*Bn(n)+(xn+(3.d0/2.d0))*Dn(n)  
c$$$ 2    continue
      do 2 n=1,Nt+1
        xn=dble(n)         
        an=(c**2)*xn*(xn+1.d0)/(sqrt(2.d0)
     &     *(2.d0*xn+1.d0)*(2.d0*xn+3.d0))
        kn=(c**2)*xn*(xn+1.d0)/(sqrt(2.d0)
     &     *(2.d0*xn+1.d0)*(2.d0*xn-1.d0))
        Dn(n)=an*(dexp((2.d0*xn+1.d0)*zetap)-dexp(2.d0*zetap))
     +        /(dcosh((2.d0*xn+1.d0)*zetap)-dcosh(2.d0*zetap)) 
        Bn(n)=kn*(dexp(-2.d0*zetap)-dexp((2.d0*xn+1.d0)*zetap))
     +        /(dcosh((2.d0*xn+1.d0)*zetap)-dcosh(2.d0*zetap))       
        dUnB(n)=(xn-5.d-1)*Bn(n)*dcosh((xn-(1.d0/2.d0))*zetap)
     &       +(xn+(3.d0/2.d0))*Dn(n)*dcosh((xn+(3.d0/2.d0))*zetap) 
        dUn(n)=(xn-5.d-1)*Bn(n)+(xn+(3.d0/2.d0))*Dn(n) 
c        write(7,*)'n=',n,'Bn=',Bn(n),'Dn=',Dn(n)
 2    continue
c      write(*,*) 'Bn',Bn
c      write(*,*) 'Dn',Dn
c      write(*,*) 'dUn',dUn 

c$$$      do 111 n=1,Nt+1
c$$$         write(*,*) 'dUn(',n,')',dUn(n)
c$$$ 111  continue   


c     Calcul des bn pour n=1,Nt
c     -------------------------
      tabbn(1)=-Bn(1)+5.d0*Dn(1)+10.d0*Bn(2)/3.d0 
c      write(*,*) 'tabbn(1)',tabbn(1)    
      do 3 n=2,Nt
       xn=dble(n)
       tabbn(n)=-(2.d0*xn-1.d0)*Bn(n)+(2.d0*xn+3.d0)*Dn(n)
     +      +(2.d0*xn*(2.d0*xn+3.d0)/(2.d0*xn+1.d0))*Bn(n+1)
     +      -(2.d0*(xn+1.d0)*(2.d0*xn-1.d0)
     +      /(2.d0*xn+1.d0))*Dn(n-1)       
c      write(*,*) 'n=',n ,'tabbn=',tabbn(n)
 3    continue
   
c     ___________________________________________________________
c     ------- Calcul des sommes alpha0 et alpha pour n=1,Nt------
      
c     Calcul de la somme sur les bm
c     -----------------------------
      sum=0.d0
      do 4 m=1,Nt-1
       xm=dble(m)       
       sum=sum+((2.d0*xm+1.d0)/(xm*(xm+1.d0)))*tabbn(m)
 4    continue  

c     On calcule les alpha0
c     ---------------------
       xNt=dble(Nt)
       alpha0=(-sum-((2.d0*xNt+1.d0)/xNt)*tabbn(Nt)) 
c       write(*,*) 'dans CalculTabs alpha0',alpha0
c     On exprime les alpha_n
c     On calcul dabord la somme sur m=1,n-1 dans lapremière boucle
c     Puis on rajoute les alpha0 et l'expression (2n+1/n)*bn
c     ------------------------------------------------------      
      alpha(2)=(3.d0/2.d0)*tabbn(1)
      
      do 5 n=3,Nt
       xn=dble(n)     
       alpha(n)=alpha(n-1)
     +  +((2.d0*(xn-1.d0)+1.d0)/((xn-1.d0)*xn))*tabbn(n-1) 
c       write(*,*) 'n=',n ,'alphan=',alpha(n)
 5    continue
      
      alpha(1)=3.d0*tabbn(1)+alpha0
c      write(*,*) 'alpha1=', alpha(1)
      do 50 n=2,Nt
       xn=dble(n)   
       alpha(n)=alpha(n)+((2.d0*xn+1.d0)/xn)*tabbn(n)+alpha0       
c       write(*,*) 'n=',n ,'alphan=',alpha(n)-alpha0
c        write(*,*) 'n=',n ,'alphan=',alpha(n)
 50   continue
      
c     ___________________________________________________________
c     ------- Calcul des S_n et S_0 -----------------------------

c     On calcul les Sn pour la surface libre et pour la surface
c     de la bulle, d'ou les tabs et tabsB

c     Expression de tabs0 et tabsB0
c     -----------------------------
      tabs0=(-alpha0-3.d0*dUn(1))
c      tabsB0=-alpha0*dcosh(zetap/2)-9.d0*dUnB(1)
      tabsB0=-3.d0*dUn(1)

c     Calcul des S_n
c     --------------
      tabs(1)=-alpha(1)-5.d0*dUn(2)+6.d0*dUn(1)
      tabsB(1)=-5.d0*dUn(2)+6.d0*dUn(1)
c      tabsB(1)=-alpha(1)*dcosh(3*zetap/2)-11.d0*dUnB(2)
c     &         +5.d0*dUnB(1)
      do 6 n=2,Nt
        xn=dble(n)
        tabs(n)=-alpha(n)-(2.d0*xn+3.d0)*dUn(n+1)
     +           +2.d0*(2.d0*xn+1.d0)*dUn(n)
     +           -(2.d0*xn-1.d0)*dUn(n-1)
c        tabsB(n)=-alpha(n)*dcosh((xn+5.d-1)*zetap)
c    &           -(2.d0*xn+9.d0)*dUnB(n+1)
c     +           +2.d0*(2.d0*xn+1.d0)*ch*dUnB(n)
c     +           -(2.d0*xn-7.d0)*dUnB(n-1)
       tabsB(n)=-(2.d0*xn+3.d0)*dUn(n+1)
     +          +2.d0*(2.d0*xn+1.d0)*dUn(n)
     +          -(2.d0*xn-1.d0)*dUn(n-1)
 6    continue


c     Fin du calcul des S_n 
c     ___________________________________________________________
c     ------- Calcul des R_n et R_0 -----------------------------

c     Ici a la difference de la sphere solide, le second membre
c     contient aussi un terme qui provient du nombre de Bond
c     et depend de lambda, c'est le Rn qui est exprime ci dessous

c     Expression de R0
c     ----------------
      R0=dexp(0.5*zetap)
      
c      write(15,*) 0, R0
c     Calcul des R_n
c     --------------
      do 7 n=1,Nt
        xn=dble(n)
        Rn(n)=(2.d0*xn+1.d0)*dexp((xn+0.5d0)*zetap)
c       write(15,*) n, Rn(n)
 7    continue   
     


c     Fin du calcul des R_n 
c      write(*,*)'On sort de CalculTabs'
      end subroutine CalculTabsB
c     -------- Fin du programme--------------------------------
c****************************************************************
