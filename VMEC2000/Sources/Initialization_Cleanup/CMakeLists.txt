
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/allocate_nunv.f
	${CMAKE_CURRENT_SOURCE_DIR}/guess_axis.f
	${CMAKE_CURRENT_SOURCE_DIR}/profil3d.f
	${CMAKE_CURRENT_SOURCE_DIR}/free_mem_ns.f
	${CMAKE_CURRENT_SOURCE_DIR}/profil1d.f
	${CMAKE_CURRENT_SOURCE_DIR}/pmass.f
	${CMAKE_CURRENT_SOURCE_DIR}/torflux.f
	${CMAKE_CURRENT_SOURCE_DIR}/free_mem_nunv.f
	${CMAKE_CURRENT_SOURCE_DIR}/piota.f.bak
	${CMAKE_CURRENT_SOURCE_DIR}/pmass.f.bak
	${CMAKE_CURRENT_SOURCE_DIR}/free_mem_funct3d.f
	${CMAKE_CURRENT_SOURCE_DIR}/fixaray.f
	${CMAKE_CURRENT_SOURCE_DIR}/allocate_funct3d.f
	${CMAKE_CURRENT_SOURCE_DIR}/load_xc_from_wout.f
	${CMAKE_CURRENT_SOURCE_DIR}/heading.f
	${CMAKE_CURRENT_SOURCE_DIR}/allocate_ns.f
	${CMAKE_CURRENT_SOURCE_DIR}/close_all_files.f
	${CMAKE_CURRENT_SOURCE_DIR}/pcurr.f
	${CMAKE_CURRENT_SOURCE_DIR}/open_output_files.f
	${CMAKE_CURRENT_SOURCE_DIR}/piota.f
	${CMAKE_CURRENT_SOURCE_DIR}/torflux_deriv.f
	${CMAKE_CURRENT_SOURCE_DIR}/vsetup.f
	${CMAKE_CURRENT_SOURCE_DIR}/free_mem_recon.f
	${CMAKE_CURRENT_SOURCE_DIR}/free_persistent_mem.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)


