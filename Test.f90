!> @author Ali Samii - 2016
!! Ali Samii - Department of ASE/EM, UT Austin
!! @brief This module is an interface between parallel Adcirc input files and ESMF library.
module ADCIRC_interpolation

    use ESMF
    use MPI

    !> \author Ali Samii - 2016
    !! \brief This object stores the data required for construction an ESMF mesh from
    !! <tt>fort.14, fort.18, partmesh.txt</tt> files.
    !!
    type meshdata
        !> vm is an ESMF_VM object.  ESMF_VM is just an ESMF virtual machine class,
        !! which we will use to get the data about the local PE and PE count.
        type(ESMF_VM)                      :: vm
        !> This array contains the node coordinates of the mesh. For
        !! example, in a 2D mesh, the \c jth coordinate of the \c nth node
        !! is stored in location <tt> 2*(n-1)+j</tt> of this array.
        real(ESMF_KIND_R8), allocatable    :: NdCoords(:)
        !> This array contains the elevation of different nodes of the mesh
        real(ESMF_KIND_R8), allocatable    :: bathymetry(:)
        !> Number of nodes present in the current PE. This is different from the
        !! number of nodes owned by this PE (cf. NumOwnedNd)
        integer(ESMF_KIND_I4)              :: NumNd
        !> Number of nodes owned by this PE. This is different from the number of
        !! nodes present in the current PE (cf. NumNd)
        integer(ESMF_KIND_I4)              :: NumOwnedNd
        !> Number of elements in the current PE. This includes ghost elements and
        !! owned elements. However, we do not bother to distinguish between owned
        !! element and present element (as we did for the nodes).
        integer(ESMF_KIND_I4)              :: NumEl
        !> Number of nodes of each element, which is simply three in 2D ADCIRC.
        integer(ESMF_KIND_I4)              :: NumND_per_El
        !> Global node numbers of the nodes which are present in the current PE.
        integer(ESMF_KIND_I4), allocatable :: NdIDs(:)
        !> Global element numbers which are present in the current PE.
        integer(ESMF_KIND_I4), allocatable :: ElIDs(:)
        !> The element connectivity array, for the present elements in the current PE.
        !! The node numbers are the local numbers of the present nodes. All the element
        !! connectivities are arranged in this one-dimensional array.
        integer(ESMF_KIND_I4), allocatable :: ElConnect(:)
        !> The number of the PE's which own each of the nodes present this PE.
        !! This number is zero-based.
        integer(ESMF_KIND_I4), allocatable :: NdOwners(:)
        !>
        integer(ESMF_KIND_I4), allocatable :: ElTypes(:)
        !> This is an array, which maps the indices of the owned nodes to the indices of the present
        !! nodes. For example, assume we are on <tt>PE = 1</tt>, and we have four nodes present, and the
        !! first and third nodes belong to <tt>PE = 0</tt>. So we have:
        !! \code
        !! NumNd = 4
        !! NumOwnedNd = 2
        !! NdOwners = (/0, 1, 0, 1/)
        !! NdIDs = (/2, 3, 5, 6/)
        !! owned_to_present = (/2, 4/)    <-- Because the first node owned by this PE is actually
        !!                                    the second node present on this PE, and so on.
        !! \endcode
        integer(ESMF_KIND_I4), allocatable :: owned_to_present_nodes(:)
    end type meshdata

contains

    !> \details As the name of this function suggests, this funciton creates a parallel
    !! ESMF_Mesh from meshdata object. This function should be called collectively by
    !! all PEs for the parallel mesh to be created. The function, extract_parallel_data_from_mesh()
    !! should be called prior to calling this function.
    !! \param the_data This the input meshdata object.
    !! \param out_esmf_mesh This is the ouput ESMF_Mesh object.
    subroutine create_parallel_esmf_mesh_from_meshdata(the_data, out_esmf_mesh)
        implicit none
        type(ESMF_Mesh), intent(out)                  :: out_esmf_mesh
        type(meshdata), intent(in)                   :: the_data
        integer, parameter                            :: dim1=2, spacedim=2, NumND_per_El=3
        integer                                       :: rc
        out_esmf_mesh=ESMF_MeshCreate(parametricDim=dim1, spatialDim=spacedim, &
            nodeIDs=the_data%NdIDs, nodeCoords=the_data%NdCoords, &
            nodeOwners=the_data%NdOwners, elementIDs=the_data%ElIDs, &
            elementTypes=the_data%ElTypes, elementConn=the_data%ElConnect, &
            rc=rc)
    end subroutine

    !> \details This function is similar to create_parallel_esmf_mesh_from_meshdata(), except that
    !! it creates a masked mesh. A masked mesh is used for example to exclude the interpolation onto
    !! some nodes, when using ESMF interpolation routines.
    !! \param in_meshdata This is the input meshdata object.
    !! \param mask_array This is an array of length NumNd (number of present nodes on this PE)
    !! which contains integer numbers. When we plan to exclude a group of nodes from interpolation,
    !! we use these mask values in the interpolation routines.
    !! \param out_masked_mesh This is the output masked ESMF_Mesh.
    subroutine create_masked_esmf_mesh_from_data(in_meshdata, mask_array, out_maked_esmf_mesh)
        implicit none
        type(ESMF_Mesh), intent(out)       :: out_maked_esmf_mesh
        type(meshdata), intent(in)        :: in_meshdata
        integer(ESMF_KIND_I4), intent(in)  :: mask_array(:)
        integer, parameter                 :: dim1=2, spacedim=2, NumND_per_El=3
        integer                            :: rc
        out_maked_esmf_mesh=ESMF_MeshCreate(parametricDim=dim1, spatialDim=spacedim, &
            nodeIDs=in_meshdata%NdIDs, nodeCoords=in_meshdata%NdCoords, &
            nodeOwners=in_meshdata%NdOwners, elementIDs=in_meshdata%ElIDs, &
            elementTypes=in_meshdata%ElTypes, elementConn=in_meshdata%ElConnect, &
            nodeMask=mask_array, rc=rc)
        !print *, "mesh with mask creation: ", rc
    end subroutine create_masked_esmf_mesh_from_data

    !> @details Using the data available in <tt> fort.14, fort.18, partmesh.txt</tt> files
    !! this function extracts the scalars and arrays required for construction of a
    !! meshdata object.
    !! After calling this fucntion, one can call create_parallel_esmf_mesh_from_meshdata()
    !! or create_masked_esmf_mesh_from_data() to create an ESMF_Mesh.
    !! @param vm This is an ESMF_VM object, which will be used to obtain the \c localPE
    !! and \c peCount of the \c MPI_Communicator.
    !! @param global_fort14_dir This is the directory path (relative to the executable
    !! or an absolute path) which contains the global \c fort.14 file (not the fort.14
    !! after decomposition).
    !! @param the_data This is the output meshdata object.
    !!
    subroutine extract_parallel_data_from_mesh(vm, global_fort14_dir, the_data)
        implicit none
        type(ESMF_VM), intent(in)             :: vm
        type(meshdata), intent(inout)        :: the_data
        character(len=*), intent(in)          :: global_fort14_dir
        character(len=6)                      :: PE_ID, garbage1
        character(len=200)                    :: fort14_filename, fort18_filename, partmesh_filename
        integer                               :: i1, j1, i_num, localPet, petCount, num_global_nodes, garbage2, garbage3
        integer, allocatable                  :: local_node_numbers(:), local_elem_numbers(:), node_owner(:)
        integer, parameter                    :: dim1=2, NumND_per_El=3

        the_data%vm = vm
        call ESMF_VMGet(vm=vm, localPet=localPet, petCount=petCount)
        write(PE_ID, "(A,I4.4)") "PE", localPet
        fort14_filename = trim(global_fort14_dir//PE_ID//"/fort.14")
        fort18_filename = trim(global_fort14_dir//PE_ID//"/fort.18")
        partmesh_filename = trim(global_fort14_dir//"/partmesh.txt")

        open(unit=14, file=fort14_filename, form='FORMATTED', status='OLD', action='READ')
        open(unit=18, file=fort18_filename, form='FORMATTED', status='OLD', action='READ')
        open(unit=100, file=partmesh_filename, form='FORMATTED', status='OLD', action='READ')

        read(unit=14, fmt=*)
        read(unit=14, fmt=*) the_data%NumEl, the_data%NumNd
        allocate(the_data%NdIDs(the_data%NumNd))
        allocate(local_node_numbers(the_data%NumNd))
        allocate(the_data%ElIDs(the_data%NumEl))
        allocate(local_elem_numbers(the_data%NumEl))
        allocate(the_data%NdCoords(dim1*the_data%NumNd))
        allocate(the_data%bathymetry(the_data%NumNd))
        allocate(the_data%ElConnect(NumND_per_El*the_data%NumEl))
        allocate(the_data%NdOwners(the_data%NumNd))
        allocate(the_data%ElTypes(the_data%NumEl))

        read(unit=18, fmt=*)
        read(unit=18, fmt=*)
        read(unit=18, fmt=*) local_elem_numbers
        the_data%ElIDs = abs(local_elem_numbers)
        read(unit=18, fmt=*) garbage1, num_global_nodes, garbage2, garbage3
        read(unit=18, fmt=*) local_node_numbers
        the_data%NumOwnedND = 0
        do i1 = 1, the_data%NumNd, 1
            if (local_node_numbers(i1) > 0) then
                the_data%NumOwnedND = the_data%NumOwnedND + 1
            end if
        end do
        the_data%NdIDs = abs(local_node_numbers)
        allocate(node_owner(num_global_nodes))
        allocate(the_data%owned_to_present_nodes(the_data%NumOwnedND))
        read(unit=100, fmt=*) node_owner

        do i1 = 1, the_data%NumNd, 1
            read(unit=14, fmt=*) local_node_numbers(i1), &
                the_data%NdCoords((i1-1)*dim1 + 1), &
                the_data%NdCoords((i1-1)*dim1 + 2), &
                the_data%bathymetry(i1)
        end do
        do i1 = 1, the_data%NumEl, 1
            read(unit=14, fmt=*) local_elem_numbers(i1), i_num, &
                the_data%ElConnect((i1-1)*NumND_per_El+1), &
                the_data%ElConnect((i1-1)*NumND_per_El+2), &
                the_data%ElConnect((i1-1)*NumND_per_El+3)
        end do

        do i1= 1, the_data%NumNd, 1
            the_data%NdOwners(i1) = node_owner(the_data%NdIDs(i1)) - 1
        end do

        j1 = 0
        do i1 = 1, the_data%NumNd, 1
            if (the_data%NdOwners(i1) == localPet) then
                j1 = j1 + 1
                the_data%owned_to_present_nodes(j1) = i1
            end if
        end do
        the_data%ElTypes = ESMF_MESHELEMTYPE_TRI

        close(14)
        close(18)
        close(100)
    end subroutine extract_parallel_data_from_mesh

    !> \details This function writes the input meshdata object to a \c vtu file.
    !! The \c vtu file is in \c XML format. This function can be used for both parallel
    !! and serial mesh writing. If one uses this function for parallel write, the
    !! processing element with \c localPE=0 should also enter this function, otherwise
    !! the \c pvtu file will not be written. This function assumes that the \c vtu file
    !! which we want to write does not exist. If we want to add fields to the files which
    !! are created before, we have to call write_node_field_to_vtu() function. If the user
    !! wants to add more data fields to the created \c vtu file, the \c last_write parameter
    !! should be passed <tt>.false.</tt> so that the program do not close the file.
    !! By closing we mean writing the last three closing lines in the XML \c vtu files.
    !! However, if this is the last time we want to write on the same \c vtu file, we
    !! have to pass \c last_write equal to <tt>.true.</tt>
    !! \param the_data This is the input data for which we create the vtu file
    !! \param vtu_filename This is the name of the vtu file
    !! \param last_write This parameter indicates if this is the last time we want to
    !! write something to this \c vtu file.
    subroutine write_meshdata_to_vtu(the_data, vtu_filename, last_write)
        implicit none
        type(meshdata), intent(in)     :: the_data
        character(len=*), intent(in)   :: vtu_filename
        integer                        :: localPet, petCount
        logical, intent(in)            :: last_write
        integer                        :: i1, indent, offset_counter, rc
        integer, parameter             :: dim1=2, spacedim=2, NumND_per_El=3, vtk_triangle=5
        indent = 0

        call ESMF_VMGet(vm=the_data%vm, localPet=localPet, petCount=petCount, rc=rc)
        if (rc .NE. ESMF_Success) then
            localPet = 0
            petCount = 1
        end if

        open(unit=1014, file=vtu_filename, form='FORMATTED', &
            status='UNKNOWN', action='WRITE')
        write(unit=1014, fmt="(A A)") '<VTKFile type="UnstructuredGrid"', &
            ' version="0.1" byte_order="BigEndian">'
        indent = indent + 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            "<UnstructuredGrid>"
        indent = indent + 2
        write(unit=1014, fmt="(A A I0 A I0 A)") repeat(" ",indent), &
            '<Piece NumberOfPoints="', the_data%NumNd, &
            '" NumberOfCells="', the_data%NumEl, '">'
        indent = indent + 2

        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '<Points>'
        indent = indent + 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
        indent = indent + 2
        do i1 = 1, the_data%NumNd, 1
            write(unit=1014, fmt="(A F0.4 ' ' F0.4 ' ' F0.4 ' ')") repeat(" ",indent), &
                the_data%NdCoords((i1-1)*dim1 + 1), &
                the_data%NdCoords((i1-1)*dim1 + 2), 0.0
        end do
        indent = indent - 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '</DataArray>'
        indent = indent - 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '</Points>'

        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '<Cells>'
        indent = indent + 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '<DataArray type="Int32" Name="connectivity" Format="ascii">'
        indent = indent + 2
        do i1 = 1, the_data%NumEl, 1
            write(unit=1014, fmt="(A I0 ' ' I0 ' ' I0 ' ')") repeat(" ",indent),&
                the_data%ElConnect((i1-1)*NumND_per_El+1)-1, &
                the_data%ElConnect((i1-1)*NumND_per_El+2)-1, &
                the_data%ElConnect((i1-1)*NumND_per_El+3)-1
        end do
        indent = indent - 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '</DataArray>'
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '<DataArray type="Int32" Name="offsets" Format="ascii">'
        indent = indent + 2
        offset_counter = 0
        do i1 = 1, the_data%NumEl, 1
            offset_counter = offset_counter + 3
            write(unit=1014, fmt="(A I0)") repeat(" ",indent), &
                offset_counter
        end do
        indent = indent - 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '</DataArray>'
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '<DataArray type="Int32" Name="types" Format="ascii">'
        indent = indent + 2
        do i1=1, the_data%NumEl, 1
            write(unit=1014, fmt="(A I2)") repeat(" ",indent), &
                vtk_triangle
        end do
        indent = indent - 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '</DataArray>'
        indent = indent - 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '</Cells>'

        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '<PointData Scalars="scalars">'
        indent = indent + 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '<DataArray type="Float32" Name="bathymetry" NumberOfComponents="1" Format="ascii">'
        indent = indent + 2
        do i1 = 1, the_data%NumNd, 1
            write(unit=1014, fmt="(A F0.4)") repeat(" ",indent), &
                the_data%bathymetry(i1)
        end do
        indent = indent - 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '</DataArray>'
        indent = indent - 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '</PointData>'

        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '<CellData Scalars="scalars">'
        indent = indent + 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '<DataArray type="Int32" Name="subdomain_id" NumberOfComponents="1" Format="ascii">'
        indent = indent + 2
        do i1 = 1, the_data%NumEl, 1
            write(unit=1014, fmt="(A I0)") repeat(" ",indent), localPet
        end do
        indent = indent - 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '</DataArray>'
        indent = indent - 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '</CellData>'
        if (last_write) then
            indent = indent - 2
            write(unit=1014, fmt="(A A)") repeat(" ",indent), &
                '</Piece>'
            indent = indent - 2
            write(unit=1014, fmt="(A A)") repeat(" ",indent), &
                '</UnstructuredGrid>'
            indent = indent - 2
            write(unit=1014, fmt="(A A)") repeat(" ",indent), &
                '</VTKFile>'
        end if
        close(1014)
    end subroutine

    !> \details This function writes the input array (\c field_array) and its name (\c field_name)
    !! to the vtu file (which should already exist and not closed). Refer to write_meshdata_to_vtu()
    !! to know more about opening vtu file and closing them. If the parameter \c last_write is true
    !! then we close this file and as such we should not write anything else on this file.
    subroutine write_node_field_to_vtu(field_array, field_name, vtu_filename, last_write)
        implicit none
        character(len=*), intent(in)   :: vtu_filename, field_name
        logical, intent(in)            :: last_write
        real(ESMF_KIND_R8), intent(in) :: field_array(:)
        integer                        :: i1, indent, num_recs
        indent = 6
        open(unit=1014, file=vtu_filename, form='FORMATTED', &
            position='APPEND', status='OLD', action='WRITE')
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '<PointData Scalars="scalars">'
        indent = indent + 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '<DataArray type="Float32" Name="'//field_name//'" NumberOfComponents="1" Format="ascii">'
        indent = indent + 2
        num_recs = size(field_array)
        do i1 = 1, num_recs, 1
            write(unit=1014, fmt="(A F0.4)") repeat(" ",indent), field_array(i1)
        end do
        indent = indent - 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '</DataArray>'
        indent = indent - 2
        write(unit=1014, fmt="(A A)") repeat(" ",indent), &
            '</PointData>'
        if (last_write) then
            indent = indent - 2
            write(unit=1014, fmt="(A A)") repeat(" ",indent), &
                '</Piece>'
            indent = indent - 2
            write(unit=1014, fmt="(A A)") repeat(" ",indent), &
                '</UnstructuredGrid>'
            indent = indent - 2
            write(unit=1014, fmt="(A A)") repeat(" ",indent), &
                '</VTKFile>'
        end if
        close(1014)
    end subroutine

    !> \details This function creates an object of type meshdata from the fort14 file given by
    !!fort14_filename. Unlike extract_parallel_data_from_mesh(), this function does not
    !! create a parallel meshdata, so it can be called by only one PE and the created meshdata
    !! object can later be used to create an ESMF_Mesh object.
    subroutine extract_global_data_from_fort14(fort14_filename, the_data)
        implicit none
        type(meshdata), intent(inout)         :: the_data
        character(len=*), intent(in)          :: fort14_filename
        integer                               :: i1, i_num
        integer, parameter                    :: dim1=2, spacedim=2, NumND_per_El=3

        open(unit=14, file=fort14_filename, form='FORMATTED', status='OLD', action='READ')
        read(unit=14, fmt=*)
        read(unit=14, fmt=*) the_data%NumEl, the_data%NumNd
        allocate(the_data%NdIDs(the_data%NumNd))
        allocate(the_data%ElIDs(the_data%NumEl))
        allocate(the_data%NdCoords(dim1*the_data%NumNd))
        allocate(the_data%bathymetry(the_data%NumNd))
        allocate(the_data%ElConnect(NumND_per_El*the_data%NumEl))
        allocate(the_data%NdOwners(the_data%NumNd))
        allocate(the_data%ElTypes(the_data%NumEl))
        do i1 = 1, the_data%NumNd, 1
            read(unit=14, fmt=*) the_data%NdIDs(i1), &
                the_data%NdCoords((i1-1)*dim1 + 1), &
                the_data%NdCoords((i1-1)*dim1 + 2), &
                the_data%bathymetry(i1)
        end do
        do i1 = 1, the_data%NumEl, 1
            read(unit=14, fmt=*) the_data%ElIDs(i1), i_num, &
                the_data%ElConnect((i1-1)*NumND_per_El+1), &
                the_data%ElConnect((i1-1)*NumND_per_El+2), &
                the_data%ElConnect((i1-1)*NumND_per_El+3)
        end do
        the_data%NdOwners = 0
        the_data%ElTypes = ESMF_MESHELEMTYPE_TRI
        close(14)
    end subroutine

    !> \details Given a local array in each PE i.e. \c fieldarray, we use MPI_Gather method to
    !! gather their elements into an array (\c out_fieldarray) in \c PE=root. For this
    !! process we use an ESMF_VM which is given to this function as an input. Since, MPI_Gather
    !! is collective this function should also be called collectively.
    subroutine gather_datafield_on_root(vm1, fieldarray, root, num_total_nodes, out_fieldarray)
        implicit none

        type(ESMF_VM), intent(in)                       :: vm1
        real(ESMF_KIND_R8), pointer, intent(in)         :: fieldarray(:)
        real(ESMF_KIND_R8), pointer, intent(out)        :: out_fieldarray(:)
        integer, intent(in)                             :: root, num_total_nodes
        integer                                         :: send_count, localPet, petCount, &
            i1, j1, k1, i_num, j_num1, rc, trash2, trash3
        integer, allocatable                            :: recv_counts(:), gather_displs(:)
        real(ESMF_KIND_R8), allocatable                 :: temp_fieldarray(:)
        character(len=6)                                :: PE_ID
        character(len=4)                                :: trash1

        call ESMF_VMGet(vm=vm1, localPet=localPet, petCount=petCount, rc=rc)
        send_count = size(fieldarray)
        if (localPet == root) then
            allocate(recv_counts(petCount))
            allocate(gather_displs(petCount))
            gather_displs(1) = 0
            recv_counts = 0
            open(unit=100, file="fine/partmesh.txt", form='FORMATTED', &
                status='OLD', action='READ')
            do i1 = 1, num_total_nodes, 1
                read(100,*) i_num
                recv_counts(i_num) = recv_counts(i_num) + 1
            end do
            do i1 = 2, petCount, 1
                gather_displs(i1) = gather_displs(i1-1) + recv_counts(i1-1)
            end do
            allocate(temp_fieldarray(num_total_nodes))
            allocate(out_fieldarray(num_total_nodes))
            close(100)
        end if

        call MPI_Gatherv(fieldarray, send_count, MPI_DOUBLE_PRECISION, &
            temp_fieldarray, recv_counts, gather_displs, MPI_DOUBLE_PRECISION, &
            0, MPI_COMM_WORLD, rc)
        print *, "gathered on root: ", rc

        if (localPet == 0) then
            do i1 = 1, petCount, 1
                write(PE_ID, "(A,I4.4)") 'PE', i1-1
                open(unit=18, file="fine/"//PE_ID//"/fort.18", form='FORMATTED', &
                    status='OLD', action='READ')
                read(unit=18, fmt=*)
                read(unit=18, fmt=*) trash1, trash2, trash3, i_num
                do j1 = 1, i_num, 1
                    read(unit=18, fmt=*)
                end do
                k1 = 0
                read(unit=18, fmt=*) trash1, trash2, trash3, i_num
                do j1 = 1, i_num, 1
                    read(unit=18, fmt=*) j_num1
                    if (j_num1 > 0) then
                        k1 = k1 + 1
                        out_fieldarray(j_num1) = temp_fieldarray(gather_displs(i1)+k1)
                    end if
                end do
            end do
        end if
    end subroutine

    !> \details This function simply deallocates the arrays created in the \c meshdata object
    !! creation steps.
    subroutine destroy_meshdata(the_data)
        implicit none
        type(meshdata), intent(inout) :: the_data
        deallocate(the_data%NdIDs)
        deallocate(the_data%ElIDs)
        deallocate(the_data%NdCoords)
        deallocate(the_data%bathymetry)
        deallocate(the_data%ElConnect)
        deallocate(the_data%NdOwners)
        deallocate(the_data%ElTypes)
    end subroutine

end module ADCIRC_interpolation

!> \mainpage
!! ## Introduction
!! We use this module to first create an object of type \ref adcirc_interpolation::meshdata
!! "meshdata" from ADCIRC input files, and then use the created meshdata object to construct an
!! ESMF_Mesh object. When we have two ESMF_Mesh objects, we can use them to interpolate
!! data between two different meshes. All of this process can be either done
!! sequentially or in parallel.
!! ### Sequential interpolation
!! For sequential interpolation we use the following path:
!!   1. Use the subroutine \ref adcirc_interpolation::extract_global_data_from_fort14()
!!      "extract_global_data_from_fort14()" to create two objects of type
!!      \ref adcirc_interpolation::meshdata "meshdata". One of these objects is used
!!      as source mesh where we want to interpolate from, and the other is used as
!!      destination mesh where we want to interpolate to.
!!   2. Create two ESMF_Mesh objects from the above \ref adcirc_interpolation::meshdata
!!      "meshdata"'s, using \ref adcirc_interpolation::create_parallel_esmf_mesh_from_meshdata
!!      "create_parallel_esmf_mesh_from_meshdata" subroutine.
!!   3. Create two \c ESMF_Field 's on the above \c ESMF_Mesh 's.
!!   4. Use ESMF library to interpolate data from the source mesh to destination mesh for
!!      those nodes that we have enough data, and extrapolate for those nodes that we do
!!      not have enough data.
!!
!! ### Parallel interpolation
!! The process is very much similar to the sequential case except in the first step we
!! use the function \ref adcirc_interpolation::extract_parallel_data_from_mesh
!! "extract_parallel_data_from_mesh" to extract mesh from the files <tt> fort.14, fort.18,
!! partmesh.txt</tt>.
!! ## Basic Usage
!! Here, we present an example of how the module adcirc_interpolation should be used.
!! Consider the following two meshes of east coast. We have decomposed each of these meshes
!! into 4 subdomain, which are shown with different colors here. The subdomain decomposition
!! is done by adcprep. The coarse mesh is used as the source mesh and the fine mesh is our
!! destination mesh. We want to interpolate the bathymetry from the source mesh to the
!! destination mesh.
!! <img src="Images/coarse_mesh_subdomains.png" width=750em />
!! <div style="text-align:center; font-size:150%;">The decomposed coarse mesh, which is used as the source mesh.</div>
!! <img src="Images/fine_mesh_subdomains.png" width=750em />
!! <div style="text-align:center; font-size:150%;">The decomposed fine mesh, which is used as the destination mesh.</div>
!! \code
!! \endcode
!!
program main

    use ESMF
    use MPI
    use ADCIRC_interpolation

    implicit none
    real(ESMF_KIND_R8), pointer                      :: mapped_field_ptr(:), unmapped_field_ptr(:), dst_maskptr(:), global_datafield(:)
    type(ESMF_VM)                                    :: vm1
    type(meshdata)                                   :: src_data, dst_data, global_src_data, global_dst_data
    type(ESMF_Mesh)                                  :: src_mesh, dst_mesh
    type(ESMF_Field)                                 :: src_field, dst_unmapped_field, dst_mapped_field, src_mask_field, dst_mask_field
    type(ESMF_RouteHandle)                           :: mapped_route_handle, unmapped_route_handle
    integer                                          :: i1, j1, k1, rc, localPet, petCount, i_num, j_num1, j_num2, send_count, trash2, trash3
    real(ESMF_KIND_R8), pointer                      :: global_field_on_root(:), sorted_global_field(:)
    integer(ESMF_KIND_I4), allocatable               :: int_dst_mask0(:), int_dst_mask1(:), &
                                                        subdomain_ID(:), recv_counts(:), gather_displs(:)
    real(ESMF_KIND_R8), allocatable                  :: src_field_array(:), mask_creator(:)
    character(len=6)                                 :: PE_ID
    character(len=4)                                 :: trash1

    character(len=:), parameter                      :: src_fort14_dir = "coarse/", dst_fort14_dir = "fine/"

    call ESMF_Initialize(vm=vm1, defaultLogFilename="test.log", &
        logKindFlag=ESMF_LOGKIND_MULTI, rc=rc)
    call ESMF_VMGet(vm=vm1, localPet=localPet, petCount=petCount, rc=rc)
    write(PE_ID, "(A,I4.4)") "PE", localPet

    call extract_parallel_data_from_mesh(vm1, src_fort14_dir, src_data)
    call create_parallel_esmf_mesh_from_meshdata(src_data, src_mesh)
    call extract_parallel_data_from_mesh(vm1, dst_fort14_dir, dst_data)
    call create_parallel_esmf_mesh_from_meshdata(dst_data, dst_mesh)
    call write_meshdata_to_vtu(src_data, PE_ID//"_src_mesh.vtu", .true.)
    call write_meshdata_to_vtu(dst_data, PE_ID//"_dst_mesh.vtu", .true.)

    if (localPet == 0) then
        call extract_global_data_from_fort14("coarse/fort.14", global_src_data)
        call write_meshdata_to_vtu(global_src_data, "coarse/global_mesh.vtu", .true.)
        call extract_global_data_from_fort14("fine/fort.14", global_dst_data)
        call write_meshdata_to_vtu(global_dst_data, "fine/global_mesh.vtu", .false.)
    end if

    allocate(mask_creator(src_data%NumOwnedNd))
    mask_creator = 1.d0
    src_mask_field = ESMF_FieldCreate(mesh=src_mesh, farray=mask_creator, &
        indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
    dst_mask_field = ESMF_FieldCreate(mesh=dst_mesh, typekind=ESMF_TYPEKIND_R8, rc=rc)
    call ESMF_FieldRegridStore(srcField=src_mask_field, dstField=dst_mask_field, &
        unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
        routeHandle=mapped_route_handle, regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
        rc=rc)
    call ESMF_FieldRegrid(srcField=src_mask_field, dstField=dst_mask_field, &
        routehandle=mapped_route_handle, rc=rc)
    call ESMF_FieldGet(dst_mask_field, farrayPtr=dst_maskptr, rc=rc)

    allocate(src_field_array(src_data%NumOwnedND))
    do i1 = 1, src_data%NumOwnedNd, 1
        src_field_array(i1) = src_data%bathymetry(src_data%owned_to_present_nodes(i1))
    end do
    src_field = ESMF_FieldCreate(mesh=src_mesh, farray= src_field_array, &
        indexflag=ESMF_INDEX_DELOCAL, rc=rc)
    dst_mapped_field = ESMF_FieldCreate(mesh=dst_mesh, typekind=ESMF_TYPEKIND_R8, rc=rc)
    dst_unmapped_field = ESMF_FieldCreate(mesh=dst_mesh, typekind=ESMF_TYPEKIND_R8, rc=rc)

    call ESMF_FieldRegridStore(srcField=src_field, dstField=dst_mapped_field, &
        unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
        routeHandle=mapped_route_handle, regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
        rc=rc)
    print *, "mapped regrid store: ", rc
    call ESMF_FieldRegridStore(srcField=src_field, dstField=dst_unmapped_field, &
        unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
        routeHandle=unmapped_route_handle, regridmethod=ESMF_REGRIDMETHOD_NEAREST_STOD, &
        rc=rc)
    print *, "unmapped regrid store: ", rc

    call ESMF_FieldRegrid(srcField=src_field, dstField=dst_mapped_field, &
        routeHandle=mapped_route_handle, rc=rc)
    print *, "mapped regriding: ", rc
    call ESMF_FieldGet(dst_mapped_field, farrayPtr=mapped_field_ptr, rc=rc)
    call ESMF_FieldRegrid(srcField=src_field, dstField=dst_unmapped_field, &
        routeHandle=unmapped_route_handle, rc=rc)
    print *, "unmapped regriding: ", rc
    call ESMF_FieldGet(dst_mapped_field, farrayPtr=mapped_field_ptr, rc=rc)
    call ESMF_FieldGet(dst_unmapped_field, farrayPtr=unmapped_field_ptr, rc=rc)

    j1 = 0
    do i1 = 1, dst_data%NumOwnedND, 1
        if (abs(dst_maskptr(i1)) < 1.d-8) then
            mapped_field_ptr(i1) = unmapped_field_ptr(i1)
        end if
    end do

    call gather_datafield_on_root(vm1, mapped_field_ptr, 0, global_dst_data%NumNd, &
        global_datafield)
    if (localPet == 0) then
        call write_node_field_to_vtu(global_datafield, "interp_bath", "fine/global_mesh.vtu", .true.)
    end if

    !deallocate(mapped_field_ptr)
    !call ESMF_FieldDestroy(src_field)

    !call ESMF_MeshDestroy(src_mesh)
    call destroy_meshdata(src_data)
    !call ESMF_MeshDestroy(dst_mesh)
    !call destroy_meshdata(dst_data)

    call ESMF_Finalize()

end program main

