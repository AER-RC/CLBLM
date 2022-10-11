module JsonModule
    implicit none
    public :: printJsonList
    public :: getJsonValue
    public :: clearJsonList
    public :: compareStrings
    public :: returnString
    public :: validateBaseKeys  ! cmb add

    interface getJsonValue
      module procedure getJsonValueBase
      module procedure getJsonValueLogical
      module procedure getJsonValueLogicalArray
      module procedure getJsonValueInteger
      module procedure getJsonValueString
      module procedure getJsonValueFloat
      module procedure getJsonValueFloatArray
      module procedure getJsonValueIntegerArray
      module procedure getJsonValueStringArray
    end interface getJsonValue

    private
    integer, public, parameter :: maxNameLength = 40
    integer, public, parameter :: maxStringLength = 100
    integer, public, parameter :: typeJson = 0
    integer, public, parameter :: typeInt  = 1
    integer, public, parameter :: typeFloat = 2
    integer, public, parameter :: typeString = 3
    integer, public, parameter :: typeLogic = 4
    integer, public, parameter :: typeNull = 5
    integer, public, parameter :: typeArray = 6

    type, public :: JsonList_t
        type(JsonElement_t)         , pointer  :: jsnObj(:) => NULL()
    end type JsonList_t

    type, public :: value_t
        integer  :: elementType
        type(JsonList_t)              :: jsnLst
        logical                       :: bulDat
        character(len=maxStringLength):: chrDat
        integer                       :: intDat
        real                          :: fltDat
        type(value_t), pointer        :: arrDat(:)
    end type value_t

    type, public :: JsonElement_t
        character(len=maxNameLength)           :: name = ''
        type(value_t)                          :: value
    end type JsonElement_t

    integer                           :: depth
    character(len=*), parameter       :: offset='                                    '
contains
  subroutine printJsonList(this)
      type(JsonList_t), intent(in) :: this
      integer   jj
      depth = 0
      print *, '>>>>>>>>>>>>>>>> JSON LIST'
      do jj = 1, size(this%jsnObj)
         call printJsonElement(this%jsnObj(jj))
      end do
      print *, '<<<<<<<<<<<<<<<< JSON LIST'
  end subroutine printJsonList

  recursive subroutine printJsonElement(this)
      type(JsonElement_t), intent(in) :: this
      write(*,'(a)', advance='no') offset(1:2*depth)//' - '// trim(this%name) // ' : '
      call printJsonValue(this%value)
  end subroutine printJsonElement

  subroutine printObject(this)
      type(JsonList_t), intent(in) :: this
      integer   jj
      do jj = 1, size(this%jsnObj)
         call printJsonElement(this%jsnObj(jj))
      end do
  end subroutine printObject

   recursive subroutine printJsonValue(dat)
     type(value_t), intent(in)  :: dat
     integer                    :: jj
     depth=depth+1
      select case(dat%elementType)
        case(typeJson)
           print *,'{'
           call printObject(dat%jsnLst)
           print *, offset(1:2*depth)//'}'
        case(typeArray)
           print *, 'Array [ '
           do jj=1, size(dat%arrDat)
              call printJsonValue(dat%arrDat(jj))
           end do
           print *, ' ]'
        case(typeString)
           print '(a)', trim(dat%chrDat)
        case(typeLogic)
           if (dat%bulDat) then
              print *,'true'
           else
              print *,'false'
           end if
        case(typeNull)
              print *,'null'
        case(typeFloat)
              print *, dat%fltDat
        case(typeInt)
              print *, dat%intDat
     end select
     depth=depth-1
   end subroutine printJsonValue

   ! CMB
   logical function findJsonBaseKey(this, key)
       integer :: i
       type(jsonlist_t), target, intent(in) :: this
       character(len=maxNameLength), intent(in) :: key
       type(jsonlist_t), pointer :: jsnlst

       findJsonBaseKey = .false.
       jsnlst => this
       jsnloop: do i = 1, size(jsnlst%jsnObj)
           if (trim(key) == trim(jsnlst%jsnobj(i)%name)) then
               findJsonBaseKey = .true.
               exit jsnloop
           endif
       enddo jsnloop
       return
   end function
   ! subroutine validateSubKeys(this, group, validSubKeys, numIsValid, errstr)
   !     character(len=maxNameLength) :: currentGroupName
   !     type(value_t), pointer :: currentGroupValue
   !     integer, intent(inout) :: numIsValid
   !
   !     currentGroupName = repeat(' ', len(currentGroupName))
   !     jsnlst => this
   !
   !     do i = 1, size(jsnlst%jsnobj)
   !         currentGroupName = trim(jsnlst%jsnobj(i)%name)
   !         currentGroupValue => jsnlst%jsnobj(i)%value
   !
   !         ! If the value is another JSON obj, we have a complex type.  This one
   !         ! should then have sub keys
   !         select case(currentGroupValue%elementType)
   !             case(typeJson)
   !                 if (associated(currentGroupValue%jsnLst%jsnObj)) then
   !                     call validateSubKeys(currentGroupValue%jsnLst, validSubKeys, isValid, errstr)
   !                     if (.not. isValid)
   !                 endif
   !             case default
   !
   !         if (currentGroupValue%elementType .eq. 0) then
   !             call validateSubKeys()
   !
   ! end subroutine
   logical function validateBaseKeys(this, validKeys, errstr)
       integer :: i, j
       type(jsonlist_t), target, intent(in) :: this
       !character(len=maxNameLength), intent(in) :: key
       character(len=*), intent(in) :: validKeys(:)
       character(len=1024), intent(out) :: errstr
       type(jsonlist_t), pointer :: jsnlst
       logical :: iPresent

       validateBaseKeys = .true.
       jsnlst => this
       jsnloop: do i = 1, size(jsnlst%jsnObj)
           iPresent = .false.
           keyloop: do j = 1, size(validKeys)
               if (trim(validKeys(j)) == trim(jsnlst%jsnobj(i)%name)) then
                   iPresent = .true.
                   exit keyloop
               endif
           enddo keyloop
           if (.not. iPresent) then
               errstr = "Group name '"//trim(jsnlst%jsnobj(i)%name)//"' is not a valid group name.  Check for typos!"
               validateBaseKeys = .false.
               return
           endif
       enddo jsnloop
   end function
   ! end CMB

   function getJsonValueBase(this, value, keyList, idx)
      type(JsonList_t), target,                   intent(in) :: this
      character(len=maxNameLength), dimension(:), intent(in) :: keyList
      integer, optional, intent(in)                          :: idx
      logical                                                :: getJsonValueBase
      type(value_t)                          , intent(inout) :: value

      integer                                                :: jj, nn, sz
      integer                                                :: kk
      type(JsonList_t), pointer :: jsnLst
      getJsonValueBase = .false.
      jsnLst => this
      do kk=1,size(keyList)
        do jj=1, size(jsnLst%jsnObj)
          if ( compareStrings(trim(keyList(kk)), trim(jsnLst%jsnObj(jj)%name)) ) then
            if (kk == size(keyList)) then
               if (present(idx)) then
                  if (jsnLst%jsnObj(jj)%value%elementType /= typeArray) then
                    sz = size(keyList)
                    print *, 'ERROR: object ', (trim(keyList(nn))//'/' , nn=1,sz-1), &
                          trim(keyList(sz)), ' is not ARRAY'
                    getJsonValueBase = .false.
                    return
                  end if
                  sz = size(jsnLst%jsnObj(jj)%value%arrDat)
                  if ( sz < idx) then
                    print *, 'ERROR: object size ', &
                          (trim(keyList(nn))// '/' , nn=1, size(keyList)-1), &
                          trim(keyList(size(keyList))), ' is ', sz
                    getJsonValueBase = .false.
                    return
                  end if

                  call copyJsonValue(jsnLst%jsnObj(jj)%value%arrDat(idx), value)
               else
                  call copyJsonValue(jsnLst%jsnObj(jj)%value, value)
               end if
               getJsonValueBase = .true.
               return
            end if
            jsnLst => jsnLst%jsnObj(jj)%value%jsnLst
            exit
          end if
        end do
      end do
   end function getJsonValueBase


   subroutine copyJsonValue(curVal, jsonVal)
      type(value_t),   intent(in)  :: curVal
      type(value_t), intent(inout) :: jsonVal
      jsonVal%elementType=curVal%elementType
      select case(curVal%elementType)
         case(typeJson)
            print *,'copying objects is not supported'
            call exit(1)
         case(typeArray)
            ! WARNING: this works, but it isn't a copy like the name implies
            ! This is needed so we can have an array of objects for
            ! solar_configuration_type.
            jsonVal%arrDat => curVal%arrDat
         case(typeString)
            jsonVal%chrDat=curVal%chrDat
         case(typeLogic)
            jsonVal%bulDat=curVal%bulDat
         case(typeFloat)
            jsonVal%fltDat=curVal%fltDat
         case(typeInt)
            jsonVal%intDat=curVal%intDat
      end select
   end subroutine copyJsonValue
!===============================================================================
! get procedures

  function getJsonValueLogical(this, logicValue,intValue, keyList, idx)
      type(JsonList_t), target,                   intent(in) :: this
      character(len=maxNameLength), dimension(:), intent(in) :: keyList
      integer, optional, intent(in)                          :: idx
      logical                                                :: getJsonValueLogical
!------>>> yma
      logical                                , intent(inout) :: logicValue
      integer                                , intent(inout) :: intValue
!------>>><<<

      type(value_t)                                          :: valueLoc

      integer                                                :: nn, sz
      integer                                                :: kk
      character(len=maxNameLength)                           :: typeName


      getJsonValueLogical = getJsonValueBase(this, valueLoc, keyList, idx)

!------>>> yma
      !--- If key not present, return value=-1
      if (.not. getJsonValueLogical) then
         !logicValue = .FALSE.
         intValue = -1
         RETURN
      endif
!------>>><<<

      if (valueLoc%elementType /= typeLogic) then
        call getTypeName(valueLoc%elementType, typeName)
        sz=size(keyList)
        print *, 'ERROR: getJsonValueLogical'
        print *, 'ERROR: object ', (trim(keyList(nn))//'/' , nn=1,sz-1), &
                trim(keyList(sz)), ' has type ', trim(typeName)
        call exit(1)
      end if

!------>>> yma
      logicValue = valueLoc%bulDat
      if (valueLoc%bulDat) then !if key == .TRUE.
         intValue =1
      else   !if key == .FALSE.
         intValue =0
      endif
!------>>><<<
  end function getJsonValueLogical

  function getJsonValueInteger(this, value, keyList, idx)
      type(JsonList_t), target,                   intent(in) :: this
      character(len=maxNameLength), dimension(:), intent(in) :: keyList
      integer, optional, intent(in)                          :: idx
      logical                                                :: getJsonValueInteger
      integer                                , intent(inout) :: value

      type(value_t)                                          :: valueLoc

      integer                                                :: nn, sz
      integer                                                :: kk
      character(len=maxNameLength)                           :: typeName


      getJsonValueInteger = getJsonValueBase(this, valueLoc, keyList, idx)
      ! CMB add: return false if couldn't find a base value
      if (.not. getJsonValueInteger) return
      ! end CMB
      if (valueLoc%elementType /= typeInt) then
        call getTypeName(valueLoc%elementType, typeName)
        sz=size(keyList)
        print *, 'ERROR: getJsonValueInteger'
        print *, 'ERROR: object ', (trim(keyList(nn))//'/' , nn=1,sz-1), &
                trim(keyList(sz)), ' has type ', trim(typeName)
        call exit(1)
      end if
      value = valueLoc%intDat
  end function getJsonValueInteger

  function getJsonValueString(this, value, keyList, idx)
      type(JsonList_t), target,                   intent(in) :: this
      character(len=maxNameLength), dimension(:), intent(in) :: keyList
      integer, optional, intent(in)                          :: idx
      logical                                                :: getJsonValueString
      character(len=*)                       , intent(inout) :: value

      type(value_t)                                          :: valueLoc

      integer                                                :: nn, sz
      integer                                                :: kk
      character(len=maxNameLength)                           :: typeName


      getJsonValueString = getJsonValueBase(this, valueLoc, keyList, idx)
      ! CMB add: return false if couldn't find a base value
      if (.not. getJsonValueString) return
      ! end CMB
      if (valueLoc%elementType /= typeString) then
        call getTypeName(valueLoc%elementType, typeName)
        sz=size(keyList)
        print *, 'ERROR: getJsonValueString'
        print *, 'ERROR: object ', (trim(keyList(nn))//'/' , nn=1,sz-1), &
                trim(keyList(sz)), ' has type ', trim(typeName)
        call exit(1)
      end if
      value = valueLoc%chrDat
  end function getJsonValueString

  function getJsonValueFloat(this, value, keyList, idx)
      type(JsonList_t), target,                   intent(in) :: this
      character(len=maxNameLength), dimension(:), intent(in) :: keyList
      integer, optional, intent(in)                          :: idx
      logical                                                :: getJsonValueFloat
      real                                   , intent(inout) :: value

      type(value_t)                                          :: valueLoc

      integer                                                :: nn, sz
      integer                                                :: kk
      character(len=maxNameLength)                           :: typeName


      getJsonValueFloat = getJsonValueBase(this, valueLoc, keyList, idx)
      ! CMB add: return false if couldn't find a base value
      if (.not. getJsonValueFloat) return
      ! end CMB
      if (valueLoc%elementType /= typeFloat) then
        call getTypeName(valueLoc%elementType, typeName)
        sz=size(keyList)
        print *, 'ERROR: getJsonValueFloat'
        print *, 'ERROR: object ', (trim(keyList(nn))//'/' , nn=1,sz-1), &
                trim(keyList(sz)), ' has type ', trim(typeName)
        call exit(1)
      end if
      value = valueLoc%fltDat
  end function getJsonValueFloat

  function getJsonValueFloatArray(this, value, keyList)
      type(JsonList_t), target,                   intent(in) :: this
      character(len=maxNameLength), dimension(:), intent(in) :: keyList
      logical                                                :: getJsonValueFloatArray
      real, dimension(:), allocatable        , intent(inout) :: value

      type(value_t)                                          :: valueLoc

      integer                                                :: kk, typeLoc
      character(len=maxNameLength)                           :: typeName
      integer                                                :: jj, nn, sz
      type(JsonList_t), pointer :: jsnLst
      getJsonValueFloatArray = .false.
      jsnLst => this
      do kk=1,size(keyList)
        do jj=1, size(jsnLst%jsnObj)
          if (trim(keyList(kk)) == trim(jsnLst%jsnObj(jj)%name)) then
            if (kk == size(keyList)) then
                if (jsnLst%jsnObj(jj)%value%elementType /= typeArray) then
                  sz = size(keyList)
                  print *, 'ERROR: object ', (trim(keyList(nn))//'/' , nn=1,sz-1), &
                        trim(keyList(sz)), ' is not ARRAY'
                  getJsonValueFloatArray = .false.
                  return
                end if
                ! check type
                typeLoc = jsnLst%jsnObj(jj)%value%arrDat(1)%elementType
                if ( typeLoc/= typeFloat) then
                  call getTypeName(typeLoc, typeName)
                  sz = size(keyList)
                  print *, 'ERROR: object ', (trim(keyList(nn))//'/' , nn=1,sz-1), &
                        trim(keyList(sz)), ' has type - ', trim(typeName)
                  getJsonValueFloatArray = .false.
                  return
                end if

                sz = size(jsnLst%jsnObj(jj)%value%arrDat)
                allocate(value(sz))
                do nn=1,sz
                  value(nn) = jsnLst%jsnObj(jj)%value%arrDat(nn)%fltDat
                end do
                getJsonValueFloatArray = .true.
                return
            end if
            jsnLst => jsnLst%jsnObj(jj)%value%jsnLst
            exit
          end if
        end do
      end do
  end function getJsonValueFloatArray

  function getJsonValueLogicalArray(this, logicValue, intValue, keyList)
      type(JsonList_t), target,                   intent(in) :: this
      character(len=maxNameLength), dimension(:), intent(in) :: keyList
      logical                                                :: getJsonValueLogicalArray
!------>>> yma
      logical, dimension(:), allocatable     , intent(inout) :: logicValue
      integer, dimension(:), allocatable     , intent(inout) :: intValue
!------>>><<<

      integer                                                :: kk, typeLoc
      character(len=maxNameLength)                           :: typeName
      integer                                                :: jj, nn, sz
      type(JsonList_t), pointer :: jsnLst

      getJsonValueLogicalArray = .false.
      jsnLst => this
      do kk=1,size(keyList)
        do jj=1, size(jsnLst%jsnObj)
          if (trim(keyList(kk)) == trim(jsnLst%jsnObj(jj)%name)) then
            if (kk == size(keyList)) then
                ! CMB: Bug fix - this is supposed to be typeArray
                !if (jsnLst%jsnObj(jj)%value%elementType /= typeLogic) then
                if (jsnLst%jsnObj(jj)%value%elementType /= typeArray) then
                  sz = size(keyList)
                  print *, 'ERROR: object ', (trim(keyList(nn))//'/' , nn=1,sz-1), &
                        trim(keyList(sz)), ' is not ARRAY'
                  getJsonValueLogicalArray = .false.
                  return
                end if
                ! check type
                typeLoc = jsnLst%jsnObj(jj)%value%arrDat(1)%elementType
                ! CMB: another bug fix - this should be type logical
                !if ( typeLoc/= typeInt) then
                if ( typeLoc/= typeLogic) then
                  call getTypeName(typeLoc, typeName)
                  sz = size(keyList)
                  print *, 'ERROR: object ', (trim(keyList(nn))//'/' , nn=1,sz-1), &
                        trim(keyList(sz)), ' has type - ', trim(typeName)
                  getJsonValueLogicalArray = .false.
                  return
                end if

                sz = size(jsnLst%jsnObj(jj)%value%arrDat)
!------>>> yma
                allocate(logicValue(sz))
                allocate(intValue(sz))
                do nn=1,sz
                  logicValue(nn) = jsnLst%jsnObj(jj)%value%arrDat(nn)%bulDat

                  if ( jsnLst%jsnObj(jj)%value%arrDat(nn)%bulDat ) then !if .TRUE. value=1 otherwise value=0
                     intValue(nn) = 1
                  else
                     intValue(nn) = 0
                  endif
                end do
!------>>><<<
                getJsonValueLogicalArray = .true.
                return
            end if
            jsnLst => jsnLst%jsnObj(jj)%value%jsnLst
            exit
          end if
        end do
      end do
  end function getJsonValueLogicalArray

  function getJsonValueIntegerArray(this, value, keyList)
      type(JsonList_t), target,                   intent(in) :: this
      character(len=maxNameLength), dimension(:), intent(in) :: keyList
      logical                                                :: getJsonValueIntegerArray
      integer, dimension(:), allocatable     , intent(inout) :: value

      integer                                                :: kk, typeLoc
      character(len=maxNameLength)                           :: typeName
      integer                                                :: jj, nn, sz
      type(JsonList_t), pointer :: jsnLst
      getJsonValueIntegerArray = .false.
      jsnLst => this
      do kk=1,size(keyList)
        do jj=1, size(jsnLst%jsnObj)
          if (trim(keyList(kk)) == trim(jsnLst%jsnObj(jj)%name)) then
            if (kk == size(keyList)) then
                if (jsnLst%jsnObj(jj)%value%elementType /= typeArray) then
                  sz = size(keyList)
                  print *, 'ERROR: object ', (trim(keyList(nn))//'/' , nn=1,sz-1), &
                        trim(keyList(sz)), ' is not ARRAY'
                  getJsonValueIntegerArray = .false.
                  return
                end if
                ! check type
                typeLoc = jsnLst%jsnObj(jj)%value%arrDat(1)%elementType
                if ( typeLoc/= typeInt) then
                  call getTypeName(typeLoc, typeName)
                  sz = size(keyList)
                  print *, 'ERROR: object ', (trim(keyList(nn))//'/' , nn=1,sz-1), &
                        trim(keyList(sz)), ' has type - ', trim(typeName)
                  getJsonValueIntegerArray = .false.
                  return
                end if

                sz = size(jsnLst%jsnObj(jj)%value%arrDat)
                allocate(value(sz))
                do nn=1,sz
                  value(nn) = jsnLst%jsnObj(jj)%value%arrDat(nn)%intDat
                end do
                getJsonValueIntegerArray = .true.
                return
            end if
            jsnLst => jsnLst%jsnObj(jj)%value%jsnLst
            exit
          end if
        end do
      end do
  end function getJsonValueIntegerArray

  function getJsonValueStringArray(this, value, keyList)
      type(JsonList_t), target,                   intent(in) :: this
      character(len=maxNameLength), dimension(:), intent(in) :: keyList
      logical                                                :: getJsonValueStringArray
      character(len=*), dimension(:), allocatable, intent(inout) :: value

      integer                                                :: kk, typeLoc
      character(len=maxNameLength)                           :: typeName
      integer                                                :: jj, nn, sz
      type(JsonList_t), pointer :: jsnLst
      getJsonValueStringArray = .false.
      jsnLst => this
      do kk=1,size(keyList)
        do jj=1, size(jsnLst%jsnObj)
          if (trim(keyList(kk)) == trim(jsnLst%jsnObj(jj)%name)) then
            if (kk == size(keyList)) then
                if (jsnLst%jsnObj(jj)%value%elementType /= typeArray) then
                  sz = size(keyList)
                  print *, 'ERROR: object ', (trim(keyList(nn))//'/' , nn=1,sz-1), &
                        trim(keyList(sz)), ' is not ARRAY'
                  getJsonValueStringArray = .false.
                  return
                end if
                ! check type
                typeLoc = jsnLst%jsnObj(jj)%value%arrDat(1)%elementType
                if ( typeLoc/= typeString) then
                  call getTypeName(typeLoc, typeName)
                  sz = size(keyList)
                  print *, 'ERROR: object ', (trim(keyList(nn))//'/' , nn=1,sz-1), &
                        trim(keyList(sz)), ' has type - ', trim(typeName)
                  getJsonValueStringArray = .false.
                  return
                end if

                sz = size(jsnLst%jsnObj(jj)%value%arrDat)
                allocate(value(sz))
                do nn=1,sz
                  value(nn) = jsnLst%jsnObj(jj)%value%arrDat(nn)%chrDat
                end do
                getJsonValueStringArray = .true.
                return
            end if
            jsnLst => jsnLst%jsnObj(jj)%value%jsnLst
            exit
          end if
        end do
      end do
  end function getJsonValueStringArray
!============================================================================
   subroutine getTypeName(curType, name)
      integer,                       intent(in)  :: curType
      character(len=maxNameLength), intent(out)  :: name
      select case(curType)
         case(typeJson)
            name = 'JSON object'
            call exit(1)
         case(typeArray)
            name = 'JSON array'
         case(typeString)
            name = 'string'
         case(typeLogic)
            name = 'logical'
         case(typeFloat)
            name = 'float'
         case(typeInt)
            name = 'integer'
      end select
   end subroutine getTypeName

!===========================================================================
  recursive subroutine clearJsonList(this)
      ! CMB: Cannot clear a const* in this context
      !CMBtype(JsonList_t), intent(in) :: this
      type(JsonList_t), intent(inout) :: this
      integer   jj
      do jj = 1, size(this%jsnObj)
        select case(this%jsnObj(jj)%value%elementType)
          case(typeJson)
             if (associated(this%jsnObj(jj)%value%jsnLst%jsnObj)) &
                      call clearJsonList(this%jsnObj(jj)%value%jsnLst)

          case(typeArray)
             deallocate(this%jsnObj(jj)%value%arrDat)
        end select
      end do
      deallocate(this%jsnObj)
  end subroutine clearJsonList

function compareStrings(strBase, strTest, caseSensitive)
   character(len=*), intent(in)  :: strBase, strTest
   logical, optional, intent(in) :: caseSensitive
   logical compareStrings
   logical caseSensitiveLoc
   character(len=maxStringLength)            :: bufBase
   character(len=maxStringLength)            :: bufTest

   if (len_trim(strBase) /= len_trim(strTest)) then
      compareStrings = .false.
      return
   endif

   if (present(caseSensitive)) then
      caseSensitiveLoc=caseSensitive
   else
      caseSensitiveLoc=.false.
   endif
   if (caseSensitiveLoc) then
      compareStrings = trim(strBase) == trim(strTest)
   else
      bufBase=strBase
      bufTest=strTest
      call toLower(bufBase)
      call toLower(bufTest)
      compareStrings = trim(bufBase) == trim(bufTest)
   endif
end function compareStrings

function returnString(strBase, caseSensitive)
   character(len=*), intent(in)  :: strBase
   logical, optional, intent(in) :: caseSensitive
   character(len=maxStringLength)  :: returnString
   logical caseSensitiveLoc

   if (present(caseSensitive)) then
      caseSensitiveLoc=caseSensitive
   else
      caseSensitiveLoc=.false.
   endif
   returnString = trim(strBase)
   if (.not. caseSensitiveLoc)  call toLower(returnString)
end function returnString

subroutine toLower(strIn)
   character(len=*), intent(inout)  :: strIn
   integer jj, kk

   do jj = 1, len_trim(strIn)
      kk = iachar(strIn(jj:jj))
      if (kk>= iachar("A") .and. kk<=iachar("Z") ) then
         strIN(jj:jj) = achar(iachar(strIn(jj:jj))+32)
      endif 
   end do
end subroutine toLower

subroutine toUpper(strIn)
   character(len=*), intent(inout)  :: strIn
   integer jj, kk

   do jj = 1, len_trim(strIn)
      kk = iachar(strIn(jj:jj))
      if (kk>= iachar("a") .and. kk<=iachar("z") ) then
         strIN(jj:jj) = achar(iachar(strIn(jj:jj))-32)
      endif 
   end do
end subroutine toUpper
end module JsonModule
