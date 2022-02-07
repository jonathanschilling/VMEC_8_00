
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/safe_open_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/control_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/vertex_list.f
	${CMAKE_CURRENT_SOURCE_DIR}/bootsj_input.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_seq.f
	${CMAKE_CURRENT_SOURCE_DIR}/vsvd0.f
	${CMAKE_CURRENT_SOURCE_DIR}/optim_params.f
	${CMAKE_CURRENT_SOURCE_DIR}/read_v3post_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/fdjac_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/de_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/vcoilpts.f
	${CMAKE_CURRENT_SOURCE_DIR}/read_boozer_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/biotsavart.f
	${CMAKE_CURRENT_SOURCE_DIR}/ga_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/saddle_surface.f
	${CMAKE_CURRENT_SOURCE_DIR}/vf_coils.f
	${CMAKE_CURRENT_SOURCE_DIR}/gade_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/bcoils_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/lmpar_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/saddle_coils.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_input.f
	${CMAKE_CURRENT_SOURCE_DIR}/mpi_params.f
	${CMAKE_CURRENT_SOURCE_DIR}/laprec.f
	${CMAKE_CURRENT_SOURCE_DIR}/coilsnamin.f
	${CMAKE_CURRENT_SOURCE_DIR}/bsc.f
	${CMAKE_CURRENT_SOURCE_DIR}/date_and_computer.f
	${CMAKE_CURRENT_SOURCE_DIR}/read_wout_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/stel_kinds.f
	${CMAKE_CURRENT_SOURCE_DIR}/fftpack.f
	${CMAKE_CURRENT_SOURCE_DIR}/modular_coils.f
	${CMAKE_CURRENT_SOURCE_DIR}/stel_constants.f
	${CMAKE_CURRENT_SOURCE_DIR}/write_array_generic.f
	${CMAKE_CURRENT_SOURCE_DIR}/vparams.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_utils.f
	${CMAKE_CURRENT_SOURCE_DIR}/mgrid_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/tf_coils.f
	${CMAKE_CURRENT_SOURCE_DIR}/system_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/bnorm_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/read_response.f
	${CMAKE_CURRENT_SOURCE_DIR}/bsc_cdf.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)



