c****************************************************************
c    Programme de constitution de P_n(x) ou P est le polynome
c     de Legendre d'ordre n et x sont argument avec |x|\leq 1
c
c    Notations:
c     n: ordrer du polynome (entier naturel positif ou nul)
c     x: argument
c     resultat: valeur de P_n(x)
c     Programmeur: A.Sellier; fevrier 2004.
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912
c     declaration des variables
      subroutine Legendre(n,x,resultat)
      integer n,k,m
      double precision x,resultat,a,xk
      dimension a(3000)
c     test
      if (n.lt.0) then
      write(*,*) 'probleme1 sur ordre du polynome de Legrendre' 
      endif
      if (n.ge.3001) then 
      write(*,*) 'probleme2 sur ordre du polynome de Legrendre' 
      endif
      if (n.eq.0) resultat=1.d0
      if (n.eq.1) resultat=x
      if (n.eq.2) resultat=(3.d0*(x**2)-1.d0)/2.d0
      if (n.ge.3) then
c      initialisation des coefficients a
       do 1 m=1,3000
        a(m)=0.d0
 1     continue
        a(1)=x
        a(2)=(3.d0*(x**2)-1.d0)/2.d0
        do 2 k=2,n-1
         xk=dble(k)
         a(k+1)=-xk*a(k-1)/(xk+1.d0)
     +          +x*(2.d0*xk+1.d0)*a(k)/(xk+1.d0)
 2      continue
        resultat=a(n)       
      endif
      end subroutine Legendre
