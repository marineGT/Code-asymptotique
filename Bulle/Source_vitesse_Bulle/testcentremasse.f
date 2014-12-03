c****************************************************************
c                Solution particuliere de l'equation 
c                   differentielle pour la surface
c                           de la bulle
c     
c
c     Appelée par la fonction gaussint1dim3 pour le calcul
c     de l'integrale cette solution pour differente valeur
c     de teta
c
c     Programmeur: M.Guemas; Fevrier 2014
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912

      function testcentremasse(teta,Nx,f,tabf,b3,e3,d3)
c     ============================================================
c     declaration des variables en entree
c     -----------------------------------
      integer ::  Nx
      double precision, intent(in)  :: teta
      double precision, dimension(Nx) :: b3, e3, d3
      double precision, dimension(0:Nx) :: f
      double precision, dimension(0:Nx) :: tabf

c     declaration des variables locales
c     ---------------------------------
      integer :: i
      double precision :: fInt

c     Declaration des variables en sorties
c     ------------------------------------
      double precision :: testcentremasse

c     Declaration des interfaces
c     ---------------------------
      interface

         function ispline(u, x, y, b, e, d, n)
           integer n
           double precision  u,x(n),y(n),b(n),e(n),d(n)
           double precision ispline
         end function ispline

      end interface
c     ================================================================== 
      
      
c     On interpole la valeur de phi correspondant au teta
c     --------------------------------------------------- 
      fInt=ispline(teta,tabf,f,b3,e3,d3,Nx) 
      
c     Calcul de la constante du CM pour une valeur de teta
c     ----------------------------------------------------
     
      testcentremasse=fInt*dsin(teta)*dcos(teta)

c      write(16,*) teta, testcentremasse

      end function testcentremasse
c     -------- Fin du programme----------------------------------
c****************************************************************

