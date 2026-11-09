"""
Running the OpenMC solver on a written model, either natively or in a container.

The container route mounts the run directory and the directory of the nuclear data at the same
absolute paths they have on the host, so the paths written into ``model.xml`` resolve the same
way inside and out, and runs as the calling user so the outputs belong to them.
"""

import glob
import hashlib
import os
import shutil
import subprocess
from pathlib import Path

from armi import runLog

from armi_openmc.options import OpenMCOptions


class OpenMCRunError(RuntimeError):
    pass


def physicalCores():
    """Number of physical cores, falling back to logical CPUs."""
    try:
        text = Path("/proc/cpuinfo").read_text()
    except OSError:
        return os.cpu_count() or 1
    cores, socket = set(), None
    for line in text.splitlines():
        if line.startswith("physical id"):
            socket = line.split(":")[1].strip()
        elif line.startswith("core id"):
            cores.add((socket, line.split(":")[1].strip()))
    return len(cores) or os.cpu_count() or 1


def resolveRunner(options: OpenMCOptions):
    """'native' or 'docker', resolving 'auto' from what is installed."""
    if options.runner != "auto":
        return options.runner
    if shutil.which(options.executable):
        return "native"
    if shutil.which("docker"):
        return "docker"
    raise OpenMCRunError(
        f"Neither the OpenMC executable `{options.executable}` nor docker was found; install one "
        "or set openmcRunner/openmcExecutable."
    )


def command(options: OpenMCOptions, workdir, args=(), program=None, mounts=(), env=None):
    """
    The command line that runs OpenMC in ``workdir``.

    Parameters
    ----------
    program : list of str, optional
        What to run instead of the ``openmc`` solver, e.g. ``["python", "script.py"]`` for a
        script that drives OpenMC through its Python API (depletion). ``args`` follow it.
    mounts : iterable of str
        More host directories the container needs, mounted read-only at their own paths.
    env : dict, optional
        More environment variables for the container.
    """
    workdir = os.path.abspath(workdir)
    threads = options.threads or physicalCores()
    runner = resolveRunner(options)
    if program is None:
        program = [options.executable if runner == "native" else "openmc", "-s", str(threads)]
    elif runner == "native" and program[0] == "python":
        program = [options.python, *program[1:]]
    if runner == "native":
        return [*program, *args], runner
    cmd = [
        "docker",
        "run",
        "--rm",
        "--name",
        containerName(workdir),
        "--user",
        f"{os.getuid()}:{os.getgid()}",
        "-e",
        "HOME=/tmp",
        "-v",
        f"{workdir}:{workdir}",
    ]
    # a multigroup model carries its data in the run directory (mgxs.h5)
    if not options.multigroup:
        dataDir = os.path.dirname(options.crossSectionsPath)
        cmd += [
            "-e",
            f"OPENMC_CROSS_SECTIONS={options.crossSectionsPath}",
            "-v",
            f"{dataDir}:{dataDir}:ro",
        ]
    for mount in mounts:
        mount = os.path.abspath(mount)
        if mount != workdir:
            cmd += ["-v", f"{mount}:{mount}:ro"]
    for key, value in (env or {}).items():
        cmd += ["-e", f"{key}={value}"]
    cmd += ["-w", workdir, options.dockerImage, *program, *args]
    return cmd, runner


def containerName(workdir):
    """A name for the container that runs ``workdir``, so an interrupted run can stop it."""
    digest = hashlib.sha1(os.path.abspath(workdir).encode()).hexdigest()[:10]
    return f"armi-openmc-{digest}-{os.getpid()}"


def run(options: OpenMCOptions, workdir, args=(), logName="openmc.out", echo=True):
    """
    Run OpenMC on the model in ``workdir`` and return the path of the last statepoint.

    The solver's screen output goes to ``workdir/openmc.out``; a failure raises with its tail.
    """
    plotting = "-p" in args or "--plot" in args
    execute(options, workdir, args=args, logName=logName, echo=echo and not plotting)
    if plotting:
        return None
    statepoint = lastStatepoint(workdir)
    if statepoint is None:
        raise OpenMCRunError(f"OpenMC wrote no statepoint in {workdir}")
    return statepoint


def execute(
    options: OpenMCOptions,
    workdir,
    args=(),
    logName="openmc.out",
    echo=True,
    program=None,
    mounts=(),
    env=None,
):
    """
    Run OpenMC, or a ``program`` that drives it, in ``workdir`` (natively or in the container).

    The screen output goes to ``workdir/<logName>`` and, when ``echo``, to ours; a failure raises
    with its tail.
    """
    workdir = os.path.abspath(workdir)
    cmd, runner = command(options, workdir, args, program=program, mounts=mounts, env=env)
    hostEnv = dict(os.environ)
    hostEnv.update(env or {})
    if not options.multigroup:
        hostEnv["OPENMC_CROSS_SECTIONS"] = options.crossSectionsPath
    plotting = "-p" in args or "--plot" in args
    runLog.info(
        f"{'Plotting the geometry with' if plotting else 'Running'} OpenMC ({runner}) in {workdir}"
    )
    runLog.debug(" ".join(cmd))
    logPath = os.path.join(workdir, logName)
    # Tee the solver's screen output: to its own log, and to ours so that whoever follows the ARMI
    # run (a terminal, the input editor's run console) sees the batches go by.
    with open(logPath, "w") as log:
        proc = subprocess.Popen(
            cmd,
            cwd=workdir,
            env=hostEnv,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        try:
            for line in proc.stdout:
                log.write(line)
                if echo:
                    print(f"[openmc] {line}", end="", flush=True)
            proc.wait()
        except BaseException:
            # An interrupted ARMI run must not leave the solver running on its own: a container
            # outlives the docker client that started it.
            if runner == "docker":
                subprocess.run(["docker", "stop", containerName(workdir)], capture_output=True)
            proc.kill()
            raise
    if proc.returncode != 0:
        tail = Path(logPath).read_text().splitlines()[-30:]
        raise OpenMCRunError(
            f"OpenMC failed in {workdir} (exit {proc.returncode}):\n" + "\n".join(tail)
        )


def lastStatepoint(workdir):
    """The newest statepoint in ``workdir``, or None."""
    statepoints = sorted(glob.glob(os.path.join(workdir, "statepoint.*.h5")), key=os.path.getmtime)
    return statepoints[-1] if statepoints else None
