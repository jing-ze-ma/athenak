"""Remove the box_convection builds of the implicit M1 tests at the end of the session
(kept when ATHENAK_M1_KEEP_BUILD is set)."""

# Modules
import os
import shutil
import pytest
import test_suite.rad_m1.m1_common as m1


@pytest.fixture(scope="session", autouse=True)
def m1_builds():
    """The binaries are built once per session (m1_common.build) and removed here."""
    yield
    if os.environ.get("ATHENAK_M1_KEEP_BUILD"):   # reuse them in the next session
        return
    for where in m1.BUILD.values():
        shutil.rmtree(where, ignore_errors=True)
