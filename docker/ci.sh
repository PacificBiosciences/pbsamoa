#!/usr/bin/env sh
set -vex

# Check if project name is provided
if [ -z "$1" ]; then
  echo "Must provide first argument with project name" 1>&2
  exit 1
fi

#options:
# To only run tests without building
: "${MESON_NINJA_BUILD:=1}"
# To build only, export MESON_RUN_TESTS=0
: "${MESON_RUN_TESTS:=1}"
# The build directory that is created and where all build products go.
: "${MESON_BUILD_DIR:=build}"

# Load correct toolchain
case "${COMPILER:-gcc}" in
  clang)
    . /app/toolchain/gcc-14/env.sh
    export CC="clang"
    export CXX="clang++"
    ;;
  *)
    . /app/toolchain/gcc-14/env.sh
    ;;
esac

# Without -fno-sanitize-recover=all UBSAN failures won't abort the program
export CFLAGS="${CFLAGS} -fno-sanitize-recover=all"
export CXXFLAGS="${CXXFLAGS} -fno-sanitize-recover=all"

# Allow to override default --werror
export MESON_WERROR="--werror"
if [ "${WERROR:-ON}" = "OFF" ]; then
  unset MESON_WERROR
fi

# Go to project directory
cd "$1"

# Link packagecache to avoid downloading packages
if [ -e subprojects ] && [ ! -e subprojects/packagecache ]; then
  ln -s /app/third-party/packagecache subprojects/packagecache
fi

# Setup project
meson setup "${MESON_BUILD_DIR}" \
  ${MESON_WERROR:-} \
  ${MESON_OPTS:-} \
  --buildtype "${BUILDTYPE:-release}" \
  --unity "${ENABLED_UNITY_BUILD:-off}" \
  --wrap-mode forcefallback \
  -Dlibdeflate:warning_level=0 \
  -Db_sanitize="${ENABLED_SANITIZERS:-none}" \
  -Dcpp_debugstl="${ENABLED_DEBUGSTL:-false}" \
  -Dtests=true

# Build project
if [ "${MESON_NINJA_BUILD}" = 1 ]; then
  ninja -C "${MESON_BUILD_DIR}" -v ${NINJA_OPTS:-}
fi

# Run tests
if [ "${MESON_RUN_TESTS}" = 1 ]; then
  meson test -C "${MESON_BUILD_DIR}" -v ${MESON_TEST_OPTS:-}
fi
