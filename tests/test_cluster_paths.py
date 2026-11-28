import os
from pathlib import Path
import shutil
import subprocess

import pytest

from smj.cluster_paths import checked_directory
from smj.errors import ConfigError, TransportError
from smj.transport import CommandResult


@pytest.mark.parametrize("root,name", [
    ("/scratch/work", "bdt"), ("~/smj-runs/sub", "bdt"),
    ("/home/user/../other/smj-runs", "bdt"), (None, "../bdt"),
    (None, "bdt--rad"), (None, "a" * 41),
])
def test_invalid_configuration_rejected_offline(root, name):
    with pytest.raises(ConfigError):
        checked_directory(None, root, name, "aims_scf")


def test_directory_resolved_by_login_shell(workdir):
    bash = shutil.which("bash")
    if os.name == "nt" and shutil.which("git"):
        bash = Path(shutil.which("git")).parents[1] / "bin/bash.exe"
    if not bash or not Path(bash).is_file():
        pytest.skip("requires Bash")

    class Shell:
        def run(self, command, *, timeout):
            p = subprocess.run([bash, "-c", command], cwd=workdir, timeout=timeout,
                               env=dict(os.environ, HOME=str(workdir.resolve())),
                               text=True, capture_output=True)
            return CommandResult(p.returncode, p.stdout, p.stderr)

    shell = Shell()
    stage = checked_directory(shell, "~/smj-runs", "TPM-B", "aims_scf")
    root = stage.rsplit("/", 2)[0]
    assert checked_directory(shell, root, "TPM-B", "aims_scf") == stage
    assert not (workdir / "smj-runs").exists()
    with pytest.raises(TransportError):
        checked_directory(shell, "/home/another/smj-runs", "TPM-B", "aims_scf")
    target = (workdir / "external").resolve()
    target.mkdir()
    run = (workdir / "smj-runs/TPM-B").resolve()
    run.parent.mkdir()
    if os.name == "nt":
        import _winapi
        _winapi.CreateJunction(str(target), str(run))
    else:
        run.symlink_to(target, target_is_directory=True)
    with pytest.raises(TransportError):
        checked_directory(shell, None, "TPM-B", "aims_scf")
