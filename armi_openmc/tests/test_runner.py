"""How OpenMC is launched."""

from armi_openmc import runner
from armi_openmc.options import OpenMCOptions
from armi_openmc.tests.conftest import crossSections, needsData


@needsData
def test_dockerContainersAreNamedAfterTheirRun(tmp_path):
    options = OpenMCOptions(runner="docker", crossSections=crossSections(), threads=2)
    cmd, kind = runner.command(options, tmp_path)
    assert kind == "docker"
    name = cmd[cmd.index("--name") + 1]
    assert name == runner.containerName(tmp_path)
    # the run directory and the data are mounted at their own paths
    assert f"{tmp_path}:{tmp_path}" in cmd
    assert cmd[-3:] == ["openmc", "-s", "2"]


def test_smearingRepresentativeBlocksAreFlagged():
    from armi import settings
    from armi.testing import mockRunLogs

    from armi_openmc.latticePhysics import warnIfSmeared

    cs = settings.Settings().modified(newSettings={"xsBlockRepresentation": "Average"})
    with mockRunLogs.BufferLog() as log:
        warnIfSmeared(cs, OpenMCOptions(blockModel="heterogeneous"))
        assert "Use Median" in log.getStdout()
    cs = cs.modified(newSettings={"xsBlockRepresentation": "Median"})
    with mockRunLogs.BufferLog() as log:
        warnIfSmeared(cs, OpenMCOptions(blockModel="heterogeneous"))
        assert "Use Median" not in log.getStdout()


def test_lastStatepointIsTheNewest(tmp_path):
    import os
    import time

    from armi_openmc import runner

    assert runner.lastStatepoint(str(tmp_path)) is None
    for n in (50, 100):
        (tmp_path / f"statepoint.{n}.h5").write_bytes(b"")
        time.sleep(0.01)
    os.utime(tmp_path / "statepoint.100.h5")
    assert runner.lastStatepoint(str(tmp_path)).endswith("statepoint.100.h5")
