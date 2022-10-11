!<f90File>**************************************************************
!
! CONTACT:
!
!   Atmospheric & Environmental Research, Inc
!   131 Hartwell Ave
!   Lexington ,MA 02421-3126 USA
!   Phone: 781.761.2288
!   E-mail: ipolonsk@aer.com
!
! COPYRIGHT NOTICE:
!
!   Copyright AER, Inc 2001-2018, All Rights Reserved
!   This software and data are covered by U.S. Patent No. 6,584,405
!   and are delivered with limited rights by AER, Inc.
!   See the file README-DATARIGHTS.txt included with this release
!   for additional details.
!
! Modified by CMB 17-Jul-2018: Enlarge buffer to 800 characters from the original 400
! Modified by Yingtao Ma 1-Aug-2018
!
!*************************************************************</f90File>

module JsonFile
   use JsonModule
   implicit none
   private
   public :: readFile

   !======================= BUFFER ====================
   !CMBcharacter(len=400), allocatable, dimension(:), target :: buffer
   !CMBcharacter(len=400), pointer, public                   :: curBuf => NULL()
   !CMBcharacter(len=400)                                    :: word
   character(len=800), allocatable, dimension(:), target :: buffer
   character(len=800), pointer, public                   :: curBuf => NULL()
   character(len=800)                                    :: word
   integer                                               :: nLine

   !======================= LIST ====================

   type, public :: listValue_t
      integer  :: valueType
      type(listObject_t), pointer   :: jsnObj(:)=> NULL() ! JSON object container
      logical                       :: bulDat
      character(len=maxStringLength):: chrDat
      integer                       :: intDat
      real                          :: fltDat
      type(listElement_t), pointer  :: jsnArr(:)=> NULL()  ! array container
   end type listValue_t

   type, public :: listElement_t
      type(listValue_t)             :: value
      type(listElement_t), pointer  :: next(:) => NULL()
   end type listElement_t

   type, public :: listObject_t
      character(len=maxNameLength)             :: name = ''
      type(listObject_t),             pointer  :: next(:)=> NULL()
      type(listValue_t)                        :: value
   end type listObject_t

! private module variables
   type(listObject_t), save, pointer           :: jsonHead(:)
   integer                                     :: depth
   logical, parameter                          :: dbg=.false.
   logical, parameter                          :: printAux=.false.

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   function readFile(fname, jsonList)
!-----------------------------------------------------------------------
      character(len=*)                :: fname
      type(JsonList_t), intent(inout) :: jsonList

      logical readFile, op
      integer, parameter  :: uid = 100
      integer             :: errStatus, jj, kk, uidNo
      !CMBcharacter(len=400), pointer  :: tempBuf
      character(len=800), pointer  :: tempBuf
      integer                      :: bufLine, bufPos

      
      inquire(file=fname, exist=readFile, opened=op, number=uidNo)
      if (.not. readFile) then
         print *, 'ERROR: ', trim(fname), ' not found'
         readFile = .false.
         call exit(1)
      end if
      if (op) close(uidNo) 

      open(uid, file=fname, form='formatted', action='read', status='old')
      nLine=1
      do
         read(uid,*,IOSTAT=errStatus)
         if (errStatus /= 0) exit
         nLine = nLine + 1
      end do
      close(uid)
      readFile = .true.
      allocate(buffer(nLine))
      
      
      open(uid, file=fname, form='formatted', action='read', status='old')
      do jj=1,nLine
         !CMBread(uid,'(A400)',IOSTAT=errStatus) buffer(jj)
         read(uid, '(A800)', IOSTAT=errStatus) buffer(jj)
      end do
      close(uid)
      if (printAux) then
         print *, '>>>>>>>>>>>>>>>>>> input file '
         print *, 'nLine', nLine
         do jj=1,nLine
            print *, trim(buffer(jj))
         end do
         print *, '<<<<<<<<<<<<<<<<<< input file '
      end if
      depth = 0
      findSymb: do bufLine=1,nLine
         bufPos=index(buffer(bufLine), '{')
         if ( bufPos /= 0 ) exit findSymb
      end do findSymb
      call incrementPos(bufLine,bufPos)
      call getObject(jsonHead,bufLine, bufPos)
      if (printAux) call printJson()
      call convertToJsonList(jsonList)
      call clearJson()
   end function readFile

! find next string between "..."
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine getName(name, bufLine, bufPos)
!-----------------------------------------------------------------------
      character(len=*), intent(inout) ::name
      integer, intent(inout) :: bufLine, bufPos
      integer                :: ps0, ln0
      ps0=bufPos
      ln0=bufLine
      call findInBuf('"', ln0, ps0)
      bufPos=ps0+1
      bufLine=ln0
      call findInBuf('"', bufLine, bufPos)
      if (bufLine==ln0) then
         curBuf => buffer(bufLine)
         name=curBuf(ps0+1:bufPos-1)
      else
         print *,'error: getName: myltiline string. fix it'
         call exit(1)
      end if
   end subroutine getName

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine findInBuf(ch, bufLine, bufPos)
!-----------------------------------------------------------------------
      character(len=1), intent(in) :: ch
      integer, intent(inout)       :: bufLine, bufPos
      integer                      :: jj, idx
      integer                      :: bufLine0, bufPos0
      bufLine0 = bufLine
      bufPos0 = bufPos
      do jj=bufLine, nLine
         curBuf => buffer(jj)
         do idx=bufPos,len_trim(curBuf)
            if (curBuf(idx:idx) == ch) then
               bufLine = jj
               bufPos = idx
               return
            end if
         end do
         bufPos=1
      end do
      print *, 'cannot find "', ch, '" in: ', bufLine0, bufPos0
      call printBufFer(bufLine0, bufPos0)
      call exit(1)
   end subroutine findInBuf

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine printBufFer(bufLine, bufPos)
!-----------------------------------------------------------------------
      integer, intent(in)       :: bufLine, bufPos
      integer                   :: jj

      curBuf => buffer(bufLine)
      print *, trim(curBuf(bufPos:len_trim(curBuf) ))

      do jj=bufLine+1,nLine
         curBuf => buffer(jj)
         print *, trim(curBuf)
      end do
   end subroutine printBufFer

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine incrementPos(bufLine, bufPos)
!-----------------------------------------------------------------------
      integer, intent(inout)    :: bufLine, bufPos
      if (len_trim(buffer(bufLine)) > bufPos) then
         bufPos=bufPos+1
      else
         bufLine=bufLine+1
         bufPos=1
      end if
   end subroutine incrementPos

! find character and stay the char position
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   function getNextCh(bufLine, bufPos)
!-----------------------------------------------------------------------
      integer, intent(inout)       :: bufLine, bufPos
      integer                      :: jj, idx
      character(len=*), parameter  :: wordStart='tfn"[{.0123456789+-'
      character(len=1)             :: getNextCh

      do bufLine=bufLine, nLine
         curBuf => buffer(bufLine)
         do bufPos=bufPos,len_trim(curBuf)
            if (index(wordStart, curBuf(bufPos:bufPos)) /= 0) then
               getNextCh=curBuf(bufPos:bufPos)
               return
            end if
         end do
         bufPos=1
      end do
      getNextCh=' '
      print *, 'failed skip space'
      call printBufFer(bufLine, bufPos)
      call exit(1)
   end function getNextCh

!-----------------------------------------------------------------------
! find word and moves to the position of the word delimeter
! "  is not word delimeter so trailing ' ' are disregarded
!-----------------------------------------------------------------------
   subroutine findWord(bufLine, bufPos)
!-----------------------------------------------------------------------
      integer, intent(inout)       :: bufLine, bufPos
      integer                      :: jj, idx
      character(len=*), parameter  :: wordDelimeter=',]}'
      integer                      :: ps

      word = ''
      ps = bufPos
      do jj=bufLine, nLine
         curBuf => buffer(jj)
         do idx=ps,len_trim(curBuf)
            if (index(wordDelimeter, curBuf(idx:idx)) /= 0) then
               word=trim(word)//curBuf(bufPos:idx-1)
               bufLine = jj
               bufPos = idx
               return
            end if
         end do
         word=trim(word)//curBuf(bufPos:len_trim(curBuf))
         ps=1
      end do
      print *, 'failed findWord'
      call printBufFer(bufLine, bufPos)
      call exit(1)
   end subroutine findWord
   !
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   recursive subroutine getValue(value, bufLine, bufPos)
!-----------------------------------------------------------------------
      type(listValue_t)              :: value
      integer, intent(inout)         :: bufLine, bufPos
      integer                        :: p0,l0
      character(len=1)               :: ch
      integer                        :: bufInt
      real                           :: bufFlt
      integer                        :: err
      ch = getNextCh(bufLine, bufPos)
      select case(ch)
         case('{')
            if (dbg) print *,'looking for Object',bufLine, bufPos
            value%valueType = typeJson
            call incrementPos(bufLine, bufPos)
            call getObject(value%jsnObj,bufLine, bufPos)
            call incrementPos(bufLine, bufPos)

         case('"')
            call incrementPos(bufLine, bufPos)
            if (dbg) print *,'looking for string' ,bufLine, bufPos
            value%valueType = typeString
            p0=bufPos
            l0=bufLine
            call findInBuf('"', bufLine, bufPos)
            if (l0==bufLine) then
               curBuf => buffer(bufLine)
               value%chrDat=curBuf(p0:bufPos-1)
            else
               print *,'multiline string. fix here'
               call exit(1)
            end if
            call incrementPos(bufLine, bufPos)

         case('t')
            call findWord(bufLine, bufPos)
            if (dbg) print *,'looking for boolian' ,bufLine, bufPos, trim(word)
            if (trim(word) == 'true') then
               value%valueType = typeLogic
               value%bulDat = .true.
            else
               print *, 'illigal entry(expected true)', trim(word), ' at ', bufLine, bufPos
               call exit(1)
            end if
            curBuf => buffer(bufLine)

         case('f')
            call findWord(bufLine, bufPos)
            if (dbg) print *,'looking for boolian' ,bufLine, bufPos, trim(word)
            if (trim(word) == 'false') then
               value%valueType = typeLogic
               value%bulDat = .false.
            else
               print *, 'illigal entry(expected false)', trim(word), ' at ', bufLine, bufPos
               call exit(1)
            end if

         case('n')
            call findWord(bufLine, bufPos)
            if (dbg) print *,'looking for null' ,bufLine, bufPos, trim(word)
            if (trim(word) == 'null') then
               value%valueType = typeNull
            else
               print *, 'illigal entry (expected null)', trim(word), ' at ', bufLine, bufPos
               call exit(1)
            end if

         case('[')   !']')
            call incrementPos(bufLine, bufPos)
            if (dbg) print *,'looking for Array' ,bufLine, bufPos
            value%valueType = typeArray
            p0=bufPos
            l0=bufLine
            call getArray(value%jsnArr,bufLine, bufPos)

            bufPos=p0
            bufLine=l0
            call findInBuf(']', bufLine, bufPos)
            curBuf => buffer(bufLine)
            value%chrDat = curBuf(p0:bufPos-1)
            if (dbg) print *,'Array : [' , trim(value%chrDat) , ' ]'
            call incrementPos(bufLine, bufPos)

         case default
            call findWord(bufLine, bufPos)
            if (dbg) print *,'looking for number' ,bufLine, bufPos, trim(word)
            if (dbg) print *, '(expected number)', trim(word), ' at ', bufLine, bufPos
            ! check if it integer
            err = index(word,'.') + index(word,'E') + index(word,'e')
            if (err /= 0) then
               value%valueType=typeFloat
               read(word,*, iostat=err) value%fltDat
            else
               value%valueType=typeInt
               read(word,*, iostat=err) value%intDat
            end if
      end select
   end subroutine getValue

!-----------------------------------------------------------------------
! subroutine finds the set of objects and saves them as a linked list
! dat - list head
! the buffer is the first character after '{'
!-----------------------------------------------------------------------
   recursive subroutine getObject(dat, bufLine, bufPos)
!-----------------------------------------------------------------------
      type(listObject_t),    pointer  :: dat(:)
      integer, intent(inout)          :: bufLine, bufPos
      character(len=1)                :: ch
      type(listObject_t),    pointer  :: curObject
      integer                         :: objectCount, cnt
      ! character(len=*), parameter     :: offset = '                                               '
      allocate(dat(1))
      depth=depth+1
      curObject=>dat(1)
      objectCount = 1
      curObject%next => NULL()
      objectLoop: do
         if (dbg) print *,'start object',bufLine, bufPos
         call getName(curObject%name, bufLine, bufPos)
         if (dbg) print *, 'processing: ', trim(curObject%name), depth, objectCount
         ! print *, 'processing: ', offset(1:4*depth), trim(curObject%name), depth, objectCount
         curBuf => buffer(bufLine)
         call findInBuf(':', bufLine, bufPos)
         call incrementPos(bufLine, bufPos)
         call getValue(curObject%value, bufLine, bufPos)

         curBuf => buffer(bufLine)
         ch = findDelimeter(bufLine, bufPos)
         if (ch == ',') then
            allocate(curObject%next(1))
            curObject => curObject%next(1)
            curObject%next => NULL()
            objectCount=objectCount+1
            call incrementPos(bufLine, bufPos)
            cycle objectLoop
         end if
         if (ch /= '}') then
            print *,'error : getObject - ', ch, bufLine, bufPos
            print *, trim(curObject%name)
            print *, 'object has to be terminated with }'
            call exit(1)
         end if
         exit objectLoop
      end do objectLoop
      depth=depth-1
   end subroutine getObject

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   recursive subroutine getArray(dat, bufLine, bufPos)
!-----------------------------------------------------------------------
      type(listElement_t),    pointer  :: dat(:)
      integer, intent(inout)         :: bufLine, bufPos
      character(len=1)               :: ch
      type(listElement_t), pointer   :: curObject
      integer                        :: objectCount
      integer                        :: lineLoc, posLoc
      allocate(dat(1))
      curObject=>dat(1)
      curObject%next => NULL()
      objectCount = 1
      curBuf => buffer(bufLine)
      arrayLoop: do
         lineLoc = bufLine
         posLoc = bufPos
         call getValue(curObject%value, bufLine, bufPos)
         curBuf => buffer(bufLine)
         ch = findDelimeter(bufLine, bufPos)
         if (ch == ',') then
            allocate(curObject%next(1))
            curObject => curObject%next(1)
            curObject%next => NULL()
            objectCount=objectCount+1
            call incrementPos(bufLine, bufPos)
            cycle arrayLoop
         end if
         if (ch /= ']') then
            print *,'error : getArray - ', ch, bufLine, bufPos
            print *, 'array has to be terminated with ]'
            call exit(1)
         end if
         exit arrayLoop
      end do arrayLoop
   end subroutine getArray

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   recursive subroutine printObject(dat)
!-----------------------------------------------------------------------
      type(listObject_t), intent(in), target  :: dat
      type(listObject_t), pointer  :: curObject
      character(len=20)            :: offset
      integer                      :: jj
      curObject => dat
      depth=depth+1
      offset=''
      offset(1:2*depth)=' '
      do
         write(*,'(a)', advance='no') offset(1:2*depth)//' - '// trim(curObject%name) // ' : '
         call printValue(curObject%value)
         if (.not. associated(curObject%next)) exit
         curObject => curObject%next(1)
      end do
      depth=depth-1
   end subroutine printObject

   ! type, public :: listValue_t
   !    integer  :: valueType
   !    type(listObject_t), pointer   :: jsnObj(:)
   !    logical                       :: bulDat
   !    character(len=maxStringLength):: chrDat
   !    integer                       :: intDat
   !    real                          :: fltDat
   !    type(listValue_t),   pointer  :: next => NULL()
   !    type(listElement_t), pointer  :: jsnArr(:)=> NULL()
   ! end type listValue_t

   ! type, public :: listElement_t
   !    type(listValue_t),   pointer  :: value(:)
   !    type(listElement_t), pointer  :: next(:) => NULL()
   ! end type listElement_t

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   recursive subroutine printValue(dat, buf)
!-----------------------------------------------------------------------
      type(listValue_t), intent(in), target  :: dat
      type(listElement_t), pointer  :: curObject
      character(len=20)             :: offset
      integer                       :: jj, off, last
      character(len=*), optional    :: buf
      select case(dat%valueType)
         case(typeJson)
            print *,'{'
            call printObject(dat%jsnObj(1))
            print *, offset(1:2*depth)//'}'
         case(typeArray)
            word=''
            last=len(word)
            off=1
            curObject => dat%jsnArr(1)
            do while(associated(curObject))
               !CMBcall printValue(curObject%value,word(off:400))
               call printValue(curObject%value,word(off:800)) 
               if (.not. associated(curObject%next)) exit
               curObject=>curObject%next(1)
               off=len_trim(word)+2
            end do
            write(*,'(a)') '[ '// trim(word)// ' ]'

         case(typeString)
            if (present(buf)) then
               write(buf,'(a)') trim(dat%chrDat)
            else
               print '(a)', trim(dat%chrDat)
            end if
         case(typeLogic)
            if (dat%bulDat) then
               if (present(buf)) then
                  write(buf,'(a)') 'true'
               else
                  print '(a)', 'true'
               end if
            else
               if (present(buf)) then
                  write(buf,'(a)') 'false'
               else
                  print '(a)', 'false'
               end if
            end if
         case(typeNull)
            if (present(buf)) then
               write(buf,'(a)') 'null'
            else
               print '(a)', 'null'
            end if
         case(typeFloat)
            if (present(buf)) then
               write(buf,*) dat%fltDat
            else
               print *, dat%fltDat
            end if
         case(typeInt)
            if (present(buf)) then
               write(buf,*) dat%intDat
            else
               print *, dat%intDat
            end if
      end select
      depth=depth-1
   end subroutine printValue

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   function findDelimeter(bufLine, bufPos)
!-----------------------------------------------------------------------
      integer, intent(inout)       :: bufLine, bufPos
      integer                      :: idx
      character(len=*), parameter  :: wordEnd=',]}:'
      integer                      :: ps
      character(len=1)             :: findDelimeter
      integer                      :: bufLine0, bufPos0

      bufLine0=bufLine
      bufPos0=bufPos

      do bufLine=bufLine,nLine
         curBuf => buffer(bufLine)
         do ps=bufPos, len_trim(curBuf)
            idx =index(wordEnd, curBuf(ps:ps))
            if (idx /= 0) then
               bufPos = ps
               findDelimeter=curBuf(ps:ps)
               return
            end if
         end do
         bufPos=1
      end do
      print *, 'failed findDelimeter',bufLine0, bufPos0
      call printBufFer(bufLine0, bufPos0)
      call exit(1)
   end function findDelimeter

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine printJson()
!-----------------------------------------------------------------------
      print *, '>>>>>>>>>>>>>>>>>> JSON OBJECT '
      if (associated(jsonHead)) then
         depth=0
         call printObject(jsonHead(1))
      end if
      print *, '<<<<<<<<<<<<<<<<<< JSON OBJECT '
   end subroutine printJson

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine convertToJsonList(jsonList)
!-----------------------------------------------------------------------
      type(JsonList_t), intent(inout) :: jsonList
      type(listObject_t), pointer  :: curObject
      integer                      :: jj
      integer                         :: cnt
      ! count object at depth
      if (associated(jsonHead)) then
         depth=0
         curObject => jsonHead(1)
         cnt=1
         do while(associated(curObject%next))
            curObject=>curObject%next(1)
            cnt=cnt+1
         end do
         allocate(jsonList%jsnObj(cnt))
         curObject => jsonHead(1)
         do jj=1,cnt
            call copyObject(curObject, jsonList%jsnObj(jj))
            if (associated(curObject%next)) curObject=>curObject%next(1)
         end do
      else
         print *,'JSON object is empty'
         return
      end if
   end subroutine convertToJsonList

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   recursive subroutine copyObject(curObject, jsonElem)
!-----------------------------------------------------------------------
      type(listObject_t),  pointer       :: curObject
      type(JsonElement_t), intent(inout) :: jsonElem
      jsonElem%name = curObject%name
      call copyValue(curObject%value, jsonElem%value)
   end subroutine copyObject

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine copyValue(curVal, jsonVal)
!-----------------------------------------------------------------------
      type(listValue_t),  intent(in)    :: curVal
      type(value_t),      intent(inout) :: jsonVal
      jsonVal%elementType=curVal%valueType
      select case(curVal%valueType)
         case(typeJson)
            call convertObject(curVal%jsnObj(1), jsonVal%jsnLst)
         case(typeArray)
            call convertArray(curVal%jsnArr(1), jsonVal%arrDat)
            jsonVal%chrDat=curVal%chrDat
         case(typeString)
            jsonVal%chrDat=curVal%chrDat
         case(typeLogic)
            jsonVal%bulDat=curVal%bulDat
         case(typeFloat)
            jsonVal%fltDat=curVal%fltDat
         case(typeInt)
            jsonVal%intDat=curVal%intDat
      end select
   end subroutine copyValue

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   recursive subroutine convertObject(dat, jsonList)
!-----------------------------------------------------------------------
      type(listObject_t), intent(in), target  :: dat
      type(JsonList_t),  intent(inout) :: jsonList

      type(listObject_t), pointer  :: curObject
      integer                      :: jj
      integer                     :: cnt
      depth=depth+1
      curObject => dat
      cnt=1
      do while(associated(curObject%next))
         cnt=cnt+1
         curObject => curObject%next(1)
      end do
      allocate(jsonList%jsnObj(cnt))
      curObject => dat
      do jj=1,cnt
         call copyObject(curObject, jsonList%jsnObj(jj))
         if (associated(curObject%next)) curObject => curObject%next(1)
      end do
      depth=depth-1
   end subroutine convertObject

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    subroutine convertArray(dat, arrDat)
!-----------------------------------------------------------------------
      type(listElement_t), intent(in), target  :: dat
      type(value_t), pointer, intent(inout)    :: arrDat(:)

      type(listElement_t), pointer  :: curObject
      integer                      :: jj
      integer                     :: cnt
      curObject => dat
      cnt=1
      do while(associated(curObject%next))
         cnt=cnt+1
         curObject => curObject%next(1)
      end do
      allocate(arrDat(cnt))
      curObject => dat
      do jj=1,cnt
         call copyValue(curObject%value, arrDat(jj))
         if (associated(curObject%next)) curObject => curObject%next(1)
      end do
   end subroutine convertArray

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine prinrBufCur(bufLine, bufPos)
!-----------------------------------------------------------------------
      integer, intent(inout)       :: bufLine, bufPos
      curBuf => buffer(bufLine)
      print *, 'prinrBufCur - |', curBuf(bufPos:bufPos),'|'
   end subroutine prinrBufCur

!-----------------------------------------------------------------------
!----------------------------------- clear memory
!-----------------------------------------------------------------------
   subroutine clearJson()
!-----------------------------------------------------------------------
      if (associated(jsonHead)) then
         call clearObject(jsonHead(1))
      end if
      deallocate(jsonHead)
   end subroutine clearJson
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   recursive subroutine clearObject(dat)
!-----------------------------------------------------------------------
      type(listObject_t), intent(in), target  :: dat
      type(listObject_t), pointer  :: curObject
      integer                      :: jj
      curObject => dat
      do
         call clearValue(curObject%value)
         if(dbg) print *, 'clearing: ', trim(curObject%name)
         if (.not. associated(curObject%next)) exit
         curObject => curObject%next(1)
      end do
   end subroutine clearObject

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   recursive subroutine clearValue(dat)
!-----------------------------------------------------------------------
      ! CMB: Cannot clear a const* in this context
      !CMBtype(listValue_t), intent(in), target  :: dat
      type(listValue_t), intent(inout), target  :: dat
      type(listElement_t), pointer  :: curObject
      type(listElement_t), pointer  :: nextObject
      character(len=20)             :: offset
      integer                       :: jj, off, last
      select case(dat%valueType)
         case(typeJson)
            call clearObject(dat%jsnObj(1))
            deallocate(dat%jsnObj)

         case(typeArray)
            curObject => dat%jsnArr(1)
            do while(associated(curObject))
               if (.not. associated(curObject%next)) exit
               nextObject => curObject%next(1)
               deallocate(curObject)
               curObject => nextObject
            end do
      end select
   end subroutine clearValue
end module JsonFile
