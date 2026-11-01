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
