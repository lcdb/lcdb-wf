import os
import pytest
import tempfile
import shutil
import inspect
from snakemake.shell import shell
from snakemake.utils import makedirs
from lcdblib.snakemake import aligners

from raw_data_fixtures import *
