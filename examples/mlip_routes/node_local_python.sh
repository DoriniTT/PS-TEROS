#!/bin/bash
# Run the Python of an environment from the node's local disk, for the GPU jobs of the MLIP route.
#
# On clusters whose shared filesystem serves many small files slowly (CephFS, NFS), importing torch and an MLIP
# package on a GPU node can take ten minutes or more per job; from local disk it takes seconds. Pack the environment
# (a self-contained conda environment, not a venv) once into a tar file that holds its folder:
#
#     tar -cf /shared/envs/my_env.tar -C /path/to/envs my_env
#
# and make this script the executable of the AiiDA code of the GPU jobs, with the tar in the code's prepend text:
#
#     export PSTEROS_ENV_TAR=/shared/envs/my_env.tar
#
# The first job on a node unpacks the tar into $PSTEROS_ENV_CACHE (default /tmp/psteros_envs); later jobs on that
# node reuse it until the tar changes (size or modification time). Jobs on the same node wait for one another while
# it unpacks. Do not replace the tar while jobs that use it are running.
set -euo pipefail

TAR=${PSTEROS_ENV_TAR:?set PSTEROS_ENV_TAR to the tar file of the environment}
CACHE=${PSTEROS_ENV_CACHE:-/tmp/psteros_envs}
NAME=$(basename "$TAR" .tar)
STAMP=$(stat -c '%s %Y' "$TAR")
mkdir -p "$CACHE"
(
  flock 9
  if [ ! -f "$CACHE/$NAME/.psteros_stamp" ] || [ "$(cat "$CACHE/$NAME/.psteros_stamp")" != "$STAMP" ]; then
    rm -rf "$CACHE/$NAME" "$CACHE/$NAME.partial"
    mkdir -p "$CACHE/$NAME.partial"
    tar -xf "$TAR" -C "$CACHE/$NAME.partial"
    test -x "$CACHE/$NAME.partial/$NAME/bin/python" || { echo "$TAR does not hold $NAME/bin/python" >&2; exit 1; }
    mv "$CACHE/$NAME.partial/$NAME" "$CACHE/$NAME"
    rmdir "$CACHE/$NAME.partial"
    echo "$STAMP" > "$CACHE/$NAME/.psteros_stamp"
  fi
) 9>"$CACHE/$NAME.lock"
exec "$CACHE/$NAME/bin/python" "$@"
