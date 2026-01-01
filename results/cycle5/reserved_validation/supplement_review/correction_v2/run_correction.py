"""Relocate the unchanged-layout review wrapper without changing its arithmetic."""
import importlib.util
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('scope_corrected_oracle',HERE/'verify_reserved.py')
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
# The v2 package is one directory deeper than the original review package.
module.WORKSPACE=HERE.parents[3]
module.REPO=module.WORKSPACE/'outputs/synonymous-rna-accessibility'
module.main()
