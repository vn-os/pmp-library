set(HEADERS MeanCurvature.h GeodesicsInHeat.h Smoothing.h)
set(SOURCES MeanCurvature.cpp GeodesicsInHeat.cpp Smoothing.cpp)

add_library(poly STATIC ${HEADERS} ${SOURCES})
target_link_libraries(poly pmp)

add_executable(polylaplace main.cpp Viewer.cpp Viewer.h)
target_link_libraries(polylaplace poly pmp_vis)

add_executable(run_tests run_tests.cpp)
target_link_libraries(run_tests poly pmp pmp_vis)
