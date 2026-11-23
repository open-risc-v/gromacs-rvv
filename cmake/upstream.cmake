# 上游接入：gromacs（cmake）→ <build>/pkg/prefix。

set(PKG_ROOT "${CMAKE_BINARY_DIR}/pkg")
if(NOT DEFINED PKG_PREFIX)
  set(PKG_PREFIX "${PKG_ROOT}/prefix")
  add_custom_target(pkg ALL
    COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/cmake/pkg-build.sh" "${SRC_DIR}" "${PKG_ROOT}"
            "$<IF:$<BOOL:${RVV_KERNELS}>,${RVV_MARCH},rv64gc>"
            "$<IF:$<BOOL:${RVV_KERNELS}>,1,0>"
    COMMENT "构建上游 gromacs（${SRC_DIR}）"
    VERBATIM)
  set(HAVE_PKG_TARGET TRUE)
endif()
