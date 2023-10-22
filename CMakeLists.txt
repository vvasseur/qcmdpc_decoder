cmake_minimum_required(VERSION 3.9)
project(qcmdpc_decoder C)
include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(qcmdpc_decoder
  src/cli.c
  src/code.c
  src/codegen.c
  src/decoder.c
  src/decoder_bp.c
  src/errorgen.c
  src/qcmdpc_decoder.c
  src/sparse_cyclic.c
  src/threshold.c
  src/xoshiro256plusplus.c)

option(AVX "Activate AVX optimization" ON)
option(PGO "Use Profile-guided optimization (set this option to GEN, then run the executable, then recompile setting this option to USE)" OFF)

if(AVX)
  target_compile_definitions(qcmdpc_decoder PUBLIC AVX=1)
endif()

foreach(option
    "PRESET_CCA"
    "PRESET_CPA"
    "INDEX"
    "BLOCK_LENGTH"
    "BLOCK_WEIGHT"
    "ERROR_WEIGHT"
    "OUROBOROS"
    "WEAK"
    "WEAK_P"
    "ERROR_FLOOR"
    "ERROR_FLOOR_P"
    "ALGO"
    "TTL_C0"
    "TTL_C1"
    "TTL_SATURATE"
    "THRESHOLD_A0"
    "THRESHOLD_A1"
    "THRESHOLD_A2"
    "THRESHOLD_A3"
    "THRESHOLD_A4"
    "BP_SCALE"
    "BP_SATURATE"
    "THRESHOLD_C0"
    "THRESHOLD_C1"
    "GRAY_SIZE")
  if(${option})
    target_compile_definitions(qcmdpc_decoder PUBLIC ${option}=${${option}})
  endif()
endforeach()

if (${CMAKE_C_COMPILER_ID} MATCHES "Clang" OR ${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast -march=native -g3")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
endif()

if (${CMAKE_C_COMPILER_ID} MATCHES "Clang")
  if(${PGO} STREQUAL "GEN")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-instr-generate=build/qcmdpc_decoder-%p.profraw")
  endif()
  if(${PGO} STREQUAL "USE")
    file(GLOB profraw_list "build/qcmdpc_decoder-*.profraw")
    execute_process(COMMAND llvm-profdata merge -output=build/qcmdpc_decoder.profdata ${profraw_list})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-instr-use=qcmdpc_decoder.profdata")
  endif()
elseif (${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
  if(${PGO} STREQUAL "GEN")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-generate=build/qcmdpc_decoder.profdata")
  endif()
  if(${PGO} STREQUAL "USE")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-use=qcmdpc_decoder.profdata -fprofile-correction")
  endif()
endif()

set_target_properties(qcmdpc_decoder
  PROPERTIES
  C_STANDARD 11
  C_STANDARD_REQUIRED YES
  C_EXTENSIONS YES
  )

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)
target_link_libraries(qcmdpc_decoder ${CMAKE_THREAD_LIBS_INIT})

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_result)
if(ipo_result)
  set_target_properties(qcmdpc_decoder PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
  target_link_libraries(qcmdpc_decoder PUBLIC ${MATH_LIBRARY})
endif()
