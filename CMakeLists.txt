include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories( $ENV{LARDATA_INC}                         )
include_directories( $ENV{LARCORE_INC}                         )
include_directories( $ENV{LARRECO_INC}                         )
include_directories( $ENV{LARSIM_INC}                          )
include_directories( $ENV{POSTGRESQL_INC}                      )

art_make(BASENAME_ONLY
         LIBRARY_NAME     EMShowerContainment
         LIB_LIBRARIES    larcore_Geometry
                          #lardata_MCBase
                          lardata_RecoBase
                          #larreco_MCComp
                          larreco_RecoAlg_ClusterRecoUtil
                          ${SIMULATIONBASE}
                          ${ART_FRAMEWORK_CORE}
                          ${ART_FRAMEWORK_PRINCIPAL}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          ${ART_FRAMEWORK_SERVICES_BASIC}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                          ${ART_PERSISTENCY_COMMON}
                          ${ART_PERSISTENCY_PROVENANCE}
                          ${ART_UTILITIES}
                          ${MF_MESSAGELOGGER}
                          ${MF_UTILITIES}
                          ${CETLIB}
                          ${ROOT_BASIC_LIB_LIST}
                          ${Boost_SYSTEM_LIBRARY}
                          ${FHICLCPP}
         MODULE_LIBRARIES larcore_Geometry
                          larcore_Geometry_Geometry_service
                          larcore_SummaryData
                          lardata_AnalysisAlg
                          lardata_AnalysisBase
                          lardata_MCBase
                          lardata_RawData
                          lardata_RecoBase
                          lardata_RecoBaseArt
                          lardata_RecoObjects
                          lardata_Utilities
                          larevt_Filters
                          larreco_ClusterFinder
                          larreco_MCComp
                          larreco_RecoAlg_ClusterRecoUtil
                          larreco_RecoAlg
                          larsim_MCCheater_BackTracker_service
                          larsim_Simulation
                          EMShowerContainment
                          ${SEEDSERVICE_SERVICE} # artextensions
                          ${SIMULATIONBASE}
                          ${ART_FRAMEWORK_CORE}
                          ${ART_FRAMEWORK_PRINCIPAL}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          ${ART_FRAMEWORK_SERVICES_BASIC}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR}
                          ${ART_PERSISTENCY_COMMON}
                          ${ART_PERSISTENCY_PROVENANCE}
                          ${ART_UTILITIES}
                          ${MF_MESSAGELOGGER}
                          ${MF_UTILITIES}
                          ${FHICLCPP}
                          ${CETLIB}
                          ${CLHEP}
                          ${ROOT_BASIC_LIB_LIST}
)

add_subdirectory(job)

install_headers()
install_fhicl()
install_source()
