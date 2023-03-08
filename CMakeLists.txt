cmake_minimum_required(VERSION 3.12)
project(QnToolsInterface CXX)
set(CMAKE_CXX_STANDARD 17)
set(PROJECT_VERSION 1.0)

find_package(ROOT REQUIRED)
find_package(QnTools REQUIRED COMPONENTS QnTools::Base QnTools::Correction)

include_directories(
        ${CMAKE_SOURCE_DIR}
        ${ROOT_INCLUDE_DIRS}
        ${QnTools_INCLUDE_DIR}/QnTools
)

message(STATUS ${ROOT_INCLUDE_DIRS})
message(STATUS ${QnTools_INCLUDE_DIR})

set(SOURCES
        src/variable_manager.cc
        src/correction_task.cc
        )
set(HEADERS
        src/variable_manager.h
        src/correction_task.h
        )

set(PCM_FILE_NAME libCorrection)
ROOT_GENERATE_DICTIONARY( G__Correction ${HEADERS}
        LINKDEF
        src/CorrectionLinkDef.h
        OPTIONS
        -I ${QnTools_INCLUDE_DIR}/QnTools
        -I ${CMAKE_SOURCE_DIR}
        )
add_library(Correction SHARED ${SOURCES} G__Correction.cxx )
target_link_libraries(Correction PUBLIC QnTools::Correction QnTools::Base ROOT::RIO ROOT::Tree ROOT::ROOTDataFrame )

add_executable( correct src/correct.cc )
target_link_libraries( correct ROOT::RIO Correction )

