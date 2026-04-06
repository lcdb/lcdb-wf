#!/usr/bin/env bash

set -euo pipefail

export DEPLOY="${DEPLOY:-/tmp/lcdb-wf-test}"
export LCDBWF_ENV="${LCDBWF_ENV:-lcdb-wf-test}"
export LCDBWF_ENV_R="${LCDBWF_ENV_R:-lcdb-wf-test-r}"
export LC_ALL="${LC_ALL:-en_US.utf8}"
export LANG="${LANG:-en_US.utf8}"
export MINIFORGE_DIR="${MINIFORGE_DIR:-$HOME/miniforge3}"

lcdbwf_repo_root() {
  git rev-parse --show-toplevel
}

lcdbwf_install_system_deps() {
  export DEBIAN_FRONTEND=noninteractive
  sudo apt-get update
  sudo apt-get install -y \
    curl \
    git \
    locales \
    locales-all \
    make \
    rsync \
    tree \
    wget \
    x11-utils
  sudo rm -rf /var/lib/apt/lists/*
  sudo localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8 || true
}

lcdbwf_init_conda() {
  if [ ! -x "${MINIFORGE_DIR}/bin/conda" ]; then
    curl -L -o /tmp/miniforge.sh \
      "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    bash /tmp/miniforge.sh -b -p "${MINIFORGE_DIR}"
  fi

  source "${MINIFORGE_DIR}/etc/profile.d/conda.sh"
  conda activate
  conda config --system --remove-key channels >/dev/null 2>&1 || true
  conda config --system --add channels bioconda
  conda config --system --add channels conda-forge
  conda config --system --set channel_priority strict
}

lcdbwf_ensure_envs() {
  lcdbwf_init_conda

  if ! conda env list | awk 'NR > 2 {print $1}' | grep -qx "${LCDBWF_ENV}"; then
    time conda env create -n "${LCDBWF_ENV}" --file env.yml
  fi

  if ! conda env list | awk 'NR > 2 {print $1}' | grep -qx "${LCDBWF_ENV_R}"; then
    time conda env create -n "${LCDBWF_ENV_R}" --file env-r.yml
  fi
}

lcdbwf_get_data() {
  local orig staging

  orig="$(lcdbwf_repo_root)"
  staging=/tmp/lcdb-wf-source

  lcdbwf_init_conda
  conda activate "${LCDBWF_ENV}"
  conda info --envs
  conda config --show

  rm -rf "${DEPLOY}" "${staging}"
  git clone "${orig}" "${staging}"
  git -C "${staging}" checkout --detach "${GITHUB_SHA:-HEAD}"

  python "${staging}/deploy.py" --flavor full --dest "${DEPLOY}"

  cp "${orig}/workflows/chipseq/run_test.sh" "${DEPLOY}/workflows/chipseq/run_test.sh"
  cp "${orig}/workflows/rnaseq/run_test.sh" "${DEPLOY}/workflows/rnaseq/run_test.sh"
  cp "${orig}/workflows/rnaseq/run_downstream_test.sh" "${DEPLOY}/workflows/rnaseq/run_downstream_test.sh"
  cp "${orig}/workflows/references/run_test.sh" "${DEPLOY}/workflows/references/run_test.sh"
  cp "${orig}/workflows/variant-calling/run_test.sh" "${DEPLOY}/workflows/variant-calling/run_test.sh"

  mkdir -p "${DEPLOY}/ci" "${DEPLOY}/test"
  cp "${orig}/test/lcdb-wf-test" "${DEPLOY}/test/lcdb-wf-test"
  cp "${orig}/test/workflow_test_params.yaml" "${DEPLOY}/test/workflow_test_params.yaml"
  cp "${orig}/ci/get-data.py" "${DEPLOY}/ci/get-data.py"
  cp "${orig}/ci/preprocessor.py" "${DEPLOY}/ci/preprocessor.py"

  (
    cd "${DEPLOY}"
    test/lcdb-wf-test data --kind=all --verbose
  )
}
