from pathlib import Path

from setuptools import setup
from setuptools.dist import Distribution


class BinaryAwareDistribution(Distribution):
    """Emit platform-specific wheels when bundled SpecHLA binaries are present."""

    def has_ext_modules(self):
        spec_hla_bin = Path(__file__).parent / "src" / "iobrpy" / "SpecHLA" / "bin"
        return spec_hla_bin.exists()


setup(distclass=BinaryAwareDistribution)
