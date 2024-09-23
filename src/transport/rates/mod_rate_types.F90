module mod_rate_types

  use precision

  implicit none

  type hydro_quantities_t
     real(kind=rk) :: rho, tem, ye, celectron, cneutron, cproton,     &
                      dneutron, dproton, delectron, dnuc, zNuc, aNuc, &
                      sumNumberDensities, sumNumberDensitiesCharge,   &
                      q_ga, effm, cneutrino
  end type

  type transport_quantities_t
     real(kind=rk) :: em
  end type

  type inelastic_scattering_opacity_t
     real(kind=rk) :: p0abr, p1abr, p2abr
     real(kind=rk) :: p0scat, p1scat
  end type

  type elastic_scattering_opacity_t
     real(kind=rk) ::  sikq, sukq
  end type elastic_scattering_opacity_t

  type absorption_opacity_t
     real(kind=rk) :: aab
  end type absorption_opacity_t

  type coupling_constants_t
     real(kind=rk) :: cc1_e(3), cc1_p(3),     &
                      cc2_e(3), cc2_p(3)
  end type


  type combined_opacities_t
     type(inelastic_scattering_opacity_t), allocatable :: inelastic_scattering(:,:,:,:)
     type(elastic_scattering_opacity_t), allocatable   :: elastic_scattering(:,:,:,:)
     type(absorption_opacity_t), allocatable           :: absorption(:,:,:,:)
     type(absorption_opacity_t), allocatable           :: kjt_absorption(:,:,:)
  end type combined_opacities_t

  interface combined_opacities_t
     module procedure contruct_combined_opacities_t
  end interface combined_opacities_t

  interface operator(+)
    module procedure add_combined_opacities_t
  end interface


  interface operator(+)
    procedure inelastic_scattering_opacity_add_scalar
    procedure elastic_scattering_opacity_add_scalar
    procedure absorption_opacity_add_scalar

!PERL-START foreach my $funcname ("inelastic_scattering", "elastic_scattering", "absorption") {
!PERL         for (my $dim = 1; $dim <=7 ; $dim++) {
    procedure @[[$funcname]]_opacity_add_@[[$dim]]d
!PERL         }
!PERL-END   }

  end interface

  contains

    function contruct_combined_opacities_t(nnmax, ie, isli, islf, t1, t2, t3, t4) result(t)

      implicit none
      integer(kind=ik), intent(in) :: nnmax, ie, isli, islf, t1, t2, t3, t4

      type(combined_opacities_t)   :: t

      type(inelastic_scattering_opacity_t), allocatable :: inelastic_scattering(:,:,:,:)
      type(elastic_scattering_opacity_t), allocatable   :: elastic_scattering(:,:,:,:)
      type(absorption_opacity_t), allocatable           :: absorption(:,:,:,:)

      type(absorption_opacity_t), allocatable           :: kjt_absorption(:,:,:)

      if (nnmax .gt. 0 .and. ie .gt. 0 .and. isli .gt. 0 .and. islf .gt. 0 &
          .and. t1 .gt. 0) then
         ! constuct type 1
         !
         ! this "should" be automatically deallocated as soon
         ! as the containing type_t t goes out of scope, right?
         !
         ! at least valgrind says it does :)
         !

         allocate(inelastic_scattering(nnmax,ie,ie,isli:islf))
         inelastic_scattering(:,:,:,:)%p0scat = 0.0_rk
         inelastic_scattering(:,:,:,:)%p1scat = 0.0_rk
         inelastic_scattering(:,:,:,:)%p0abr  = 0.0_rk
         inelastic_scattering(:,:,:,:)%p1abr  = 0.0_rk
         inelastic_scattering(:,:,:,:)%p2abr  = 0.0_rk

      endif


      if (nnmax .gt. 0 .and. ie .gt. 0 .and. isli .gt. 0 .and. islf .gt. 0 &
          .and. t2 .gt. 0) then
         ! constuct type 1
         !
         ! this "should" be automatically deallocated as soon
         ! as the containing type_t t goes out of scope, right?
         !
         ! at least valgrind says it does :)
         !

         allocate(elastic_scattering(nnmax,ie,isli:islf,3))
         elastic_scattering(:,:,:,:)%sikq = 0.0_rk
         elastic_scattering(:,:,:,:)%sukq = 0.0_rk

      endif

      if (nnmax .gt. 0 .and. ie .gt. 0 .and. isli .gt. 0 .and. islf .gt. 0 &
           .and. t3 .gt. 0) then
         allocate(absorption(nnmax,ie,isli:islf,5))
         absorption(:,:,:,:)%aab = 0.0_rk
      endif

      if (nnmax .gt. 0 .and. ie .gt. 0 .and. isli .gt. 0 .and. islf .gt. 0 &
          .and. t4 .gt. 0) then
         allocate(kjt_absorption(nnmax,ie,isli:islf))
         kjt_absorption(:,:,:)%aab = 0.0_rk
      endif

      t = combined_opacities_t(inelastic_scattering, elastic_scattering, &
                               absorption, kjt_absorption)

    end function contruct_combined_opacities_t

    function add_combined_opacities_t(type1, type2) result(added)
      use precision
      use abort
      implicit none

      type(combined_opacities_t)             :: added
      type(combined_opacities_t), intent(in) :: type1, type2

      integer(kind=ik) :: inelastic_shape(4), elastic_shape(4), &
                          absorption_shape(4), kjt_absorption_shape(3)
      integer(kind=ik) :: t1, t2, t3, t4
      integer(kind=ik) :: nnmax, ie, isi, isf

      logical          :: specialAdd


      specialAdd = .false.

      t1 = 0
      t2 = 0
      t3 = 0
      t4 = 0

      nnmax = 0
      ie    = 0
      isi   = 0
      isf   = 0

      ! check inelastic component
      inelastic_shape(:) = 0
      if (allocated(type1%inelastic_scattering)) then
         t1 = 1
         inelastic_shape(:) = shape(type1%inelastic_scattering)

         nnmax = size(type1%inelastic_scattering,dim=1)
         ie    = size(type1%inelastic_scattering,dim=2)
         isi   = lbound(type1%inelastic_scattering,dim=4)
         isf   = ubound(type1%inelastic_scattering,dim=4)
      endif
      if (allocated(type2%inelastic_scattering)) then
         t1 = 1
         if (all(inelastic_shape(:) == 0)) then
            inelastic_shape(:) = shape(type2%inelastic_scattering)

            nnmax = size(type1%inelastic_scattering,dim=1)
            ie    = size(type1%inelastic_scattering,dim=2)
            isi   = lbound(type1%inelastic_scattering,dim=4)
            isf   = ubound(type1%inelastic_scattering,dim=4)

         else if (any(inelastic_shape(:) /= shape(type2%inelastic_scattering))) then

            ! special case for nunu reactions


            if  ((      size(type1%inelastic_scattering,dim=1)       &
                 .eq. size(type2%inelastic_scattering,dim=1) .and. &
                      size(type1%inelastic_scattering,dim=2)       &
                 .eq. size(type2%inelastic_scattering,dim=2) .and. &
                      size(type1%inelastic_scattering,dim=3)       &
                 .eq. size(type2%inelastic_scattering,dim=3) .and. &
                      lbound(type2%inelastic_scattering,dim=4) .eq. 3)) then

               nnmax = size(type1%inelastic_scattering,dim=1)
               ie    = size(type1%inelastic_scattering,dim=2)
               isi   = lbound(type1%inelastic_scattering,dim=4)
               isf   = ubound(type1%inelastic_scattering,dim=4)


               if (isf .gt. 3) then
                  raise_abort("no more neutrinos_implemted yet")
               endif

               specialAdd = .true.

            else
               print *,"Cannot add two types with incompatible dimensions in the inelastic compontent togehter"
               print *,"size(type1%inelastic_scattering,dim=1) = ", size(type1%inelastic_scattering,dim=1)
               print *,"size(type2%inelastic_scattering,dim=1) = ", size(type2%inelastic_scattering,dim=1)
               print *,"size(type1%inelastic_scattering,dim=2) = ", size(type1%inelastic_scattering,dim=2)
               print *,"size(type2%inelastic_scattering,dim=2) = ", size(type2%inelastic_scattering,dim=2)
               print *,"size(type1%inelastic_scattering,dim=3) = ", size(type1%inelastic_scattering,dim=3)
               print *,"size(type2%inelastic_scattering,dim=3) = ", size(type2%inelastic_scattering,dim=3)
               print *,"lbound(type1%inelastic_scattering,dim=4) = ",lbound(type1%inelastic_scattering,dim=4)
               print *,"lbound(type2%inelastic_scattering,dim=4) = ",lbound(type2%inelastic_scattering,dim=4)
               raise_abort("add_combined_opacities_t")
            endif
         endif
      endif

      ! check elastic component
      elastic_shape(:) = 0
      if (allocated(type1%elastic_scattering)) then
         t2 = 1
         elastic_shape(:) = shape(type1%elastic_scattering)


         nnmax = size(type1%elastic_scattering,dim=1)
         ie    = size(type1%elastic_scattering,dim=2)
         isi   = lbound(type1%elastic_scattering,dim=3)
         isf   = ubound(type1%elastic_scattering,dim=3)
      endif
      if (allocated(type2%elastic_scattering)) then
         t2 = 1
         if (all(elastic_shape(:) == 0)) then
            elastic_shape(:) = shape(type2%elastic_scattering)

            nnmax = size(type1%elastic_scattering,dim=1)
            ie    = size(type1%elastic_scattering,dim=2)
            isi   = lbound(type1%elastic_scattering,dim=3)
            isf   = ubound(type1%elastic_scattering,dim=3)

         else if (any(elastic_shape(:) /= shape(type2%elastic_scattering))) then
            print *,"Cannot add two types with incompatible dimensions in the elastic compontent togehter"
               print *,"size(type1%elastic_scattering,dim=1) = ", size(type1%elastic_scattering,dim=1)
               print *,"size(type2%elastic_scattering,dim=1) = ", size(type2%elastic_scattering,dim=1)
               print *,"size(type1%elastic_scattering,dim=2) = ", size(type1%elastic_scattering,dim=2)
               print *,"size(type2%elastic_scattering,dim=2) = ", size(type2%elastic_scattering,dim=2)
               print *,"size(type1%elastic_scattering,dim=3) = ", size(type1%elastic_scattering,dim=3)
               print *,"size(type2%elastic_scattering,dim=3) = ", size(type2%elastic_scattering,dim=3)
               print *,"size(type1%elastic_scattering,dim=4) = ", size(type1%elastic_scattering,dim=4)
               print *,"size(type2%elastic_scattering,dim=4) = ", size(type2%elastic_scattering,dim=4)
            raise_abort("add_combined_opacities_t")
         endif
      endif
!!$
!!$      ! check absorption component
      absorption_shape(:) = 0
      if (allocated(type1%absorption)) then
         t3 = 1
         absorption_shape(:) = shape(type1%absorption)

         nnmax = size(type1%absorption,dim=1)
         ie    = size(type1%absorption,dim=2)
         isi   = lbound(type1%absorption,dim=3)
         isf   = ubound(type1%absorption,dim=3)
      endif
      if (allocated(type2%absorption)) then
         t3 = 1
         if (all(absorption_shape(:) == 0)) then
            absorption_shape(:) = shape(type2%absorption)

            nnmax = size(type1%absorption,dim=1)
            ie    = size(type1%absorption,dim=2)
            isi   = lbound(type1%absorption,dim=3)
            isf   = ubound(type1%absorption,dim=3)
         else if (any(absorption_shape(:) /= shape(type2%absorption))) then
            print *,"Cannot add two types with incompatible dimensions in the absoprtion compontent togehter"
            raise_abort("add_combined_opacities_t")
         endif
      endif
!!$
      ! check kjt_absorption component
      kjt_absorption_shape(:) = 0
      if (allocated(type1%kjt_absorption)) then
         t4 = 1
         kjt_absorption_shape(:) = shape(type1%kjt_absorption)

         nnmax = size(type1%kjt_absorption,dim=1)
         ie    = size(type1%kjt_absorption,dim=2)
         isi   = lbound(type1%kjt_absorption,dim=3)
         isf   = ubound(type1%kjt_absorption,dim=3)
      endif
      if (allocated(type2%kjt_absorption)) then
         t4 = 1
         if (all(kjt_absorption_shape(:) == 0)) then
            kjt_absorption_shape(:) = shape(type2%kjt_absorption)

            nnmax = size(type1%kjt_absorption,dim=1)
            ie    = size(type1%kjt_absorption,dim=2)
            isi   = lbound(type1%kjt_absorption,dim=3)
            isf   = ubound(type1%kjt_absorption,dim=3)
         else if (any(kjt_absorption_shape(:) /= shape(type2%kjt_absorption))) then
            print *,"Cannot add two types with incompatible dimensions in the kjt_absoprtion compontent togehter"
            raise_abort("add_combined_opacities_t")
         endif
      endif

      added = combined_opacities_t(nnmax,ie,isi,isf,t1,t2,t3,t4)

      if (.not.(specialAdd)) then

         if (allocated(type1%inelastic_scattering)) then
            added%inelastic_scattering(:,:,:,:)%p0scat =         &
                 added%inelastic_scattering(:,:,:,:)%p0scat    + &
                 type1%inelastic_scattering(:,:,:,:)%p0scat

            added%inelastic_scattering(:,:,:,:)%p1scat =         &
                 added%inelastic_scattering(:,:,:,:)%p1scat    + &
                 type1%inelastic_scattering(:,:,:,:)%p1scat

            added%inelastic_scattering(:,:,:,:)%p0abr =         &
                 added%inelastic_scattering(:,:,:,:)%p0abr    + &
                 type1%inelastic_scattering(:,:,:,:)%p0abr

            added%inelastic_scattering(:,:,:,:)%p1abr =         &
                 added%inelastic_scattering(:,:,:,:)%p1abr    + &
                 type1%inelastic_scattering(:,:,:,:)%p1abr

            added%inelastic_scattering(:,:,:,:)%p2abr =         &
                 added%inelastic_scattering(:,:,:,:)%p2abr    + &
                 type1%inelastic_scattering(:,:,:,:)%p2abr

         endif

         if (allocated(type2%inelastic_scattering)) then
            added%inelastic_scattering(:,:,:,:)%p0scat =         &
                 added%inelastic_scattering(:,:,:,:)%p0scat    + &
                 type2%inelastic_scattering(:,:,:,:)%p0scat

            added%inelastic_scattering(:,:,:,:)%p1scat =         &
                 added%inelastic_scattering(:,:,:,:)%p1scat    + &
                 type2%inelastic_scattering(:,:,:,:)%p1scat

            added%inelastic_scattering(:,:,:,:)%p0abr =         &
                 added%inelastic_scattering(:,:,:,:)%p0abr    + &
                 type2%inelastic_scattering(:,:,:,:)%p0abr

            added%inelastic_scattering(:,:,:,:)%p1abr =         &
                 added%inelastic_scattering(:,:,:,:)%p1abr    + &
                 type2%inelastic_scattering(:,:,:,:)%p1abr

            added%inelastic_scattering(:,:,:,:)%p2abr =         &
                 added%inelastic_scattering(:,:,:,:)%p2abr    + &
                 type2%inelastic_scattering(:,:,:,:)%p2abr

         endif


      else ! specialAdd

 !        print *,"performing special Add"

         if (allocated(type1%inelastic_scattering)) then
            added%inelastic_scattering(:,:,:,:)%p0scat =         &
                 added%inelastic_scattering(:,:,:,:)%p0scat    + &
                 type1%inelastic_scattering(:,:,:,:)%p0scat

            added%inelastic_scattering(:,:,:,:)%p1scat =         &
                 added%inelastic_scattering(:,:,:,:)%p1scat    + &
                 type1%inelastic_scattering(:,:,:,:)%p1scat

            added%inelastic_scattering(:,:,:,:)%p0abr =         &
                 added%inelastic_scattering(:,:,:,:)%p0abr    + &
                 type1%inelastic_scattering(:,:,:,:)%p0abr

            added%inelastic_scattering(:,:,:,:)%p1abr =         &
                 added%inelastic_scattering(:,:,:,:)%p1abr    + &
                 type1%inelastic_scattering(:,:,:,:)%p1abr

            added%inelastic_scattering(:,:,:,:)%p2abr =         &
                 added%inelastic_scattering(:,:,:,:)%p2abr    + &
                 type1%inelastic_scattering(:,:,:,:)%p2abr

         endif


         if (allocated(type2%inelastic_scattering)) then
            added%inelastic_scattering(:,:,:,3:3)%p0scat =         &
                 added%inelastic_scattering(:,:,:,3:3)%p0scat    + &
                 type2%inelastic_scattering(:,:,:,3:3)%p0scat

            added%inelastic_scattering(:,:,:,3:3)%p1scat =         &
                 added%inelastic_scattering(:,:,:,3:3)%p1scat    + &
                 type2%inelastic_scattering(:,:,:,3:3)%p1scat

            added%inelastic_scattering(:,:,:,3:3)%p0abr =         &
                 added%inelastic_scattering(:,:,:,3:3)%p0abr    + &
                 type2%inelastic_scattering(:,:,:,3:3)%p0abr

            added%inelastic_scattering(:,:,:,3:3)%p1abr =         &
                 added%inelastic_scattering(:,:,:,3:3)%p1abr    + &
                 type2%inelastic_scattering(:,:,:,3:3)%p1abr

            added%inelastic_scattering(:,:,:,3:3)%p2abr =         &
                 added%inelastic_scattering(:,:,:,3:3)%p2abr    + &
                 type2%inelastic_scattering(:,:,:,3:3)%p2abr

         endif

      endif ! specialAdd


      if (allocated(type1%elastic_scattering)) then
         added%elastic_scattering(:,:,:,:)%sikq =         &
              added%elastic_scattering(:,:,:,:)%sikq    + &
              type1%elastic_scattering(:,:,:,:)%sikq

         added%elastic_scattering(:,:,:,:)%sukq =         &
              added%elastic_scattering(:,:,:,:)%sukq    + &
              type1%elastic_scattering(:,:,:,:)%sukq
      endif

      if (allocated(type2%elastic_scattering)) then
         added%elastic_scattering(:,:,:,:)%sikq =         &
              added%elastic_scattering(:,:,:,:)%sikq    + &
              type2%elastic_scattering(:,:,:,:)%sikq

         added%elastic_scattering(:,:,:,:)%sukq =         &
              added%elastic_scattering(:,:,:,:)%sukq    + &
              type2%elastic_scattering(:,:,:,:)%sukq
      endif

      if (allocated(type1%absorption)) then
         added%absorption(:,:,:,:)%aab =         &
              added%absorption(:,:,:,:)%aab    + &
              type1%absorption(:,:,:,:)%aab
      endif

      if (allocated(type2%absorption)) then
         added%absorption(:,:,:,:)%aab =         &
              added%absorption(:,:,:,:)%aab    + &
              type2%absorption(:,:,:,:)%aab
      endif

      if (allocated(type1%kjt_absorption)) then
         added%kjt_absorption(:,:,:)%aab =         &
              added%kjt_absorption(:,:,:)%aab    + &
              type1%kjt_absorption(:,:,:)%aab
      endif

      if (allocated(type2%kjt_absorption)) then
         added%kjt_absorption(:,:,:)%aab =         &
              added%kjt_absorption(:,:,:)%aab    + &
              type2%kjt_absorption(:,:,:)%aab
      endif

    end function add_combined_opacities_t



!PERL-START my %components = (
!PERL         "inelastic_scattering" => ["p0abr", "p1abr", "p2abr", "p0scat", "p1scat"],
!PERL         "elastic_scattering" => ["sikq", "sukq"],
!PERL         "absorption" => ["aab"]
!PERL       );
!PERL       foreach my $funcname ("inelastic_scattering", "elastic_scattering", "absorption") {

    pure function @[[$funcname]]_opacity_add_scalar(opacity1, opacity2) result(add)
      type(@[[$funcname]]_opacity_t), intent(in) :: opacity1, opacity2
      type(@[[$funcname]]_opacity_t) :: add

!PERL       foreach my $component (@{$components{$funcname}}) {
      add%@[[$component]] = opacity1%@[[$component]] + opacity2%@[[$component]]
!PERL       }
    end function

!PERL       for (my $dim = 1; $dim <= 7; $dim++) {
!PERL         my $colons = join(",", ((":") x $dim));
!PERL         my $opac1shape = join(", &\n    ", map { "size(opacity1, dim=$_)" } 1..$dim);
    function @[[$funcname]]_opacity_add_@[[$dim]]d(opacity1, opacity2) result(add)
      type(@[[$funcname]]_opacity_t), intent(in) :: &
        opacity1(@[[$colons]]), opacity2(@[[$opac1shape]])
      type(@[[$funcname]]_opacity_t) :: &
        add(@[[$opac1shape]])

!PERL           foreach my $component (@{$components{$funcname}}) {
      add(@[[$colons]])%@[[$component]] = &
          opacity1(@[[$opac1shape]])%@[[$component]] &
        + opacity2(@[[$opac1shape]])%@[[$component]]

!PERL           }

    end function

!PERL        }
!PERL-END }

end module mod_rate_types

#ifdef PROGRAM_operator_test
program operator_test
!make NO_CONFIG=1
  use precision
  use mod_rate_types
  use abort

  implicit none
  type(inelastic_scattering_opacity_t) :: a(1), b(1), c(1)
  type(absorption_opacity_t) :: aa(1:3,1), ab(3:5,1), ac(1:5,1)

  ! test normal case
  a(1)%p0abr = 0.0_rk
  a(1)%p1abr = 1.0_rk

  b(1)%p0abr = 2.0_rk
  b(1)%p1abr = 3.0_rk

  c = a + b

  abort_if(c(1)%p0abr /= 2.0_rk)
  abort_if(c(1)%p1abr /= 4.0_rk)

  ! test non-trivial overlap
  aa(:,:)%aab = 1.0_rk
  ab(:,:)%aab = 2.0_rk
  ac = aa + ab
  abort_if(ac(1,1)%aab /= 1.0_rk)
  abort_if(ac(2,1)%aab /= 1.0_rk)
  abort_if(ac(3,1)%aab /= 3.0_rk)
  abort_if(ac(4,1)%aab /= 2.0_rk)
  abort_if(ac(5,1)%aab /= 2.0_rk)
end program
#endif
